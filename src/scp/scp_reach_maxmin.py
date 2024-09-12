from __future__ import division
import stormpy
import stormpy.core
import stormpy.logic
import stormpy.pars
import re
import stormpy.examples
import stormpy.examples.files
import time
from gurobipy import *
from pipenv.patched.notpip.cmdoptions import timeout
#import scp.interval_parser as interval_parser

from util import interval_parser


from collections import defaultdict
import collections
import stormpy.info

import pycarl
import pycarl.core



import stormpy.pomdp

import stormpy._config as config

if not config.storm_with_pars:
    print("Support parameters is missing. Try building storm-pars.")
    raise AttributeError

import stormpy.pars
from pycarl.formula import FormulaType, Relation

if stormpy.info.storm_ratfunc_use_cln():
    import pycarl.cln.formula
else:
    import pycarl.gmp.formula
pycarl.clear_pools()

def tree():
    return collections.defaultdict(tree)

# Utiliy function to create dictionary
def multi_dict(K, type):
    if K == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: multi_dict(K - 1, type))


def input_files():
    path = "collision_partial_obs_2d_upd_hobs_20_small.prism"
    interval_path="collision_partial_obs_2d_upd_hobs_20_big.intervals"
    formula_str = "P=? [F \"goal\"]"
    threshold = 0.75
    return path, interval_path, formula_str, threshold



class QcqpOptions():
    def __init__(self, mu, maxiter, graph_epsilon, silent, timeout):
        self.mu = mu
        self.maxiter = maxiter
        self.graph_epsilon = graph_epsilon
        self.silent = silent
        self.timeout = timeout


class QcqpResult():
    def __init__(self, value_at_initial, parameter_values):
        self.value_at_initial = value_at_initial
        self.parameter_values = parameter_values




#class for incremental encoding
class QcqpSolver_affine_simple_fun():
    def __init__(self):
        self.solver_timer = 0.0
        self.encoding_timer = 0.0
        self.robust_encoding_timer = 0.0
        self.model_check_timer =0.0

        self.constraint_timer = 0.0
        self._constants_floats = dict()
        self._encoding = None

        self.iterations = 0
        self._pVars = None
        self._tau = None
        self._paramVars = None
        self._dual_upper = None
        self._dual_lower = None
        self._mu=None
        self.solver_params = None
        self.solver_output= []

    def _float_repr(self, constant_val):
        """
        Returns the float representation for a constant value
        :param constant_val:
        :return:
        """
        if constant_val.is_one():
            return 1.0
        elif constant_val.is_minus_one():
            return -1.0

        v = self._constants_floats.get(constant_val, float(constant_val))
        self._constants_floats[constant_val] = v
        return v

    def _create_encoding(self,model):
        """
        Creates the encoding for Gurobi
        :param model: stormpy model for uPOMDP
        :return:
        """
        numstate = model.nr_states

        self._encoding = Model("qcp")
        self._encoding.setParam('OutputFlag', not self._options.silent)

        # Initializing some arrays for state, parameter and penalty variables, and their values at previous iterations
        # Initializing gurobi variables for parameters,lb=lowerbound, ub=upperbound
        self._pVars = [self._encoding.addVar(lb=0, ub=1.0) for _ in range(numstate)]
        self._tau = [self._encoding.addVar(lb=0) for _ in range(numstate)]

        # Initializing policy variables
        self._paramVars = dict([[x.id, self._encoding.addVar(lb=self._options.graph_epsilon, ub=1 - self._options.graph_epsilon)] for x in self._parameters])
        self._dual_upper = dict()
        self._dual_lower = dict()

        # Updates the model for gurobi
        numcons = 0

        # Initializes the dual variables for each interval
        for item in self._intervals:
            stateval = int(item.state)
            succval = int(item.successor)
            self._dual_upper[(stateval, succval)] = self._encoding.addVar(lb=0)
            self._dual_lower[(stateval, succval)] = self._encoding.addVar(lb=0)
        self._encoding.update()

        #Gurobi parameters
        self._encoding.Params.OutputFlag = 0
        self._encoding.Params.Presolve = 2
        self._encoding.Params.Method = 2
        self._encoding.Params.Crossover = 0
        self._encoding.Params.CrossoverBasis = 0
        self._encoding.Params.NumericFocus = 3
        self._encoding.Params.BarHomogeneous = 1
        self._encoding.Params.ScaleFlag = 3
        self._encoding.Params.FeasibilityTol = 1e-6
        self._encoding.Params.OptimalityTol = 1e-6
        self._encoding.Params.BarConvTol = 1e-6
        self._encoding.update()

    def _model_constraints(self,model,i):
        """
        Creates the model constraints
        :param model: stormpy model for uPOMDP
        :param i: iteration number
        :return:
        """
        # start encoding
        for state in model.states:
            #if int(state.id) % 10000 == 0:
            #    print("Encoding through states",state.id)
            start2 = time.time()
            #for first iteration, build all states
            if i==1:
                #if a state with action, it is a non-robust constraint
                if int(state.id) not in self._robuststates:
                    start2 = time.time()
                    # non robust constraint, there are only 1 actions in our encoding:
                    for action in state.actions:
                        cons = 0
                        #for each transition in a state
                        for transition in action.transitions:
                            #gather the transitioning state
                            succ = int(transition.column)
                            # Value of transition
                            transition_value = transition.value()
                            # denominator of the value of transition
                            den = transition_value.denominator.constant_part()

                            # If the transition value is not constant
                            if not transition_value.is_constant():
                                num = transition_value.numerator.polynomial()

                                # Iterates over terms in numerators
                                for t in num:
                                    coeff = self._float_repr(t.coeff / den)

                                    # If the transition term is a constant
                                    if t.is_constant():
                                        # Add just value of transition is the successor state is prob1 to the constraints
                                        if self._prob1A.get(succ):
                                            cons = cons + coeff
                                        # Else add transitionvalue*p_successor to the constraint
                                        else:
                                            cons = cons + coeff * self._pVars[succ]

                                    # If the transition term is not a constant, which means
                                    # the probability of taking an action
                                    else:
                                        # the degree should be 1 for simple pomdps
                                        if t.tdeg > 1:
                                            raise RuntimeError("We expect the term to be a single variable")
                                        # interval parameter should not be in a state with actions for simple pomdps
                                        if t.monomial[0][0] in self._interval_parameters:
                                            raise RuntimeError("Interval parameter in non-robust constraint")
                                        #gather the parameter that corresponds to an action
                                        param_id = t.monomial[0][0].id

                                        # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1
                                        if self._prob1A.get(succ):
                                            cons = cons + coeff * self._paramVars[param_id]
                                        # Otherwise, add the convexified terms instead of bilinear constraints
                                        else:
                                            cons = cons + coeff * self._pinit[succ] * (self._paramVars[param_id])
                                            cons = cons + coeff * self._paraminit[param_id] * (self._pVars[succ] - self._pinit[succ])

                            # If the transition value not constant
                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if self._prob1A.get(succ):
                                    cons = cons + self._float_repr(constant_value)
                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + self._float_repr(constant_value) * self._pVars[succ]

                        try:
                            # additional step to make the model more numerically stable by removing
                            # small enough constants, they should be accompained by the penalty variables
                            if abs(cons.getConstant()) < 1e-6:
                                cons = cons - cons.getConstant()
                        except:
                            pass
                        self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] <= cons + self._tau[state.id]))
                    end2 = time.time()
                    self.robust_encoding_timer += (end2 - start2)
                # for other iterations, build states with multiple actions after each iteration
                else:

                    start2 = time.time()
                    # non robust constraint
                    for action in state.actions:
                        cons = 0
                        #for each transition in a state

                        for transition in action.transitions:
                            #gather the transitioning state
                            succ = int(transition.column)
                            # Value of transition
                            transition_value = transition.value()
                            # Denominator of transition
                            den = transition_value.denominator.constant_part()

                            # If the transition value is not constant
                            if not transition_value.is_constant():
                                num = transition_value.numerator.polynomial()

                                # Iterates over terms in numerators
                                for t in num:
                                    coeff = self._float_repr(t.coeff / den)

                                    # If the transition term is a constant
                                    if t.is_constant():
                                        # Add just value of transition is the successor state is prob1 to the constraints
                                        if self._prob1A.get(succ):
                                            cons = cons + coeff
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            cons = cons + coeff * self._pVars[succ]

                                    # If the transition term is not a constant
                                    else:
                                        #check if the degree is 1, it should be 1 for simple pomdps
                                        if t.tdeg == 1:

                                            #gather the uncertain interval in the state
                                            if t.monomial[0][0] in self._interval_parameters:
                                                interval_name = t.monomial[0][0].name
                                                interval_value = -1
                                                for item in self._items:
                                                    if item.name == interval_name:
                                                        interval_value = 0
                                                        break
                                                if interval_value == -1:
                                                    # print("interval name: " + str(interval_name))
                                                    raise RuntimeError("Interval not part of the instantiation")
                                                stateval = int(state.id)
                                                succval = int(succ)
                                                # add robust constraint with dualization
                                                self._encoding.addConstr(self._dual_upper[(stateval, succval)] - self._dual_lower[
                                                    (stateval, succval)] == coeff * self._pVars[succ])
                                                # add terms to the state constraints
                                                cons = cons +self._dual_upper[(stateval, succval)] * item.upperbound - \
                                                       self._dual_lower[(stateval, succval)] * item.lowerbound
                                            #if the parameter is an policy parameter:
                                            elif t.monomial[0][0] in self._parameters:
                                                # gather parameter id
                                                param_id = t.monomial[0][0].id
                                                # Add the value of policy parameter if the probability of successor state is 1
                                                if self._prob1A.get(succ):
                                                    cons = cons + coeff * paramVars[param_id]
                                                # Add nothing successor state is prob0
                                                else:
                                                    # add convexified term to the constraints
                                                    cons = cons + coeff * self._pinit[succ] * (self._paramVars[param_id])
                                                    cons = cons + coeff * self._paraminit[param_id] * (self._pVars[succ] - self._pinit[succ])
                                        #for simple pomdps, this degree should be 1
                                        elif t.tdeg > 1:
                                            raise RuntimeError("the supplied POMDP is not simple")
                            # If the value of transition is constant
                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if prob1A.get(succ):
                                    cons = cons + self._float_repr(constant_value)

                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + self._float_repr(constant_value) * self._pVars[succ]
                        try:
                            # additional step to make the model more numerically stable by removing
                            # small enough constants, they should be accompained by the penalty variables
                            if abs(cons.getConstant()) < 1e-6:
                                cons = cons - cons.getConstant()
                        except:
                            pass
                        self._encoding.addConstr(self._pVars[state.id] <= cons + self._tau[state.id])
                    end2 = time.time()
                    self.robust_encoding_timer += (end2 - start2)
            #for other iterations, just build the states with actions
            else:

                if int(state.id) not in self._robuststates:
                    start2 = time.time()

                    # non robust constraint:
                    for action in state.actions:
                        cons = 0
                        for transition in action.transitions:

                            succ = int(transition.column)
                            # Value of transition
                            transition_value = transition.value()
                            den = transition_value.denominator.constant_part()
                            # denom1 = 1 / self._float_repr(den)

                            # If the transition value is not constant
                            if not transition_value.is_constant():
                                num = transition_value.numerator.polynomial()

                                # Iterates over terms in numerators
                                for t in num:
                                    coeff = self._float_repr(t.coeff / den)

                                    # If the transition term is a constant
                                    if t.is_constant():
                                        # Add just value of transition is the successor state is prob1 to the constraints
                                        if self._prob1A.get(succ):
                                            cons = cons + coeff
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            cons = cons + coeff * self._pVars[succ]

                                    # If the transition term is not a constant
                                    else:
                                        # the degree should be 1 for simple pomdps
                                        if t.tdeg > 1:
                                            raise RuntimeError("We expect the term to be a single variable")
                                        # interval parameter should not be in a state with actions for simple pomdps
                                        if t.monomial[0][0] in self._interval_parameters:
                                            raise RuntimeError("Interval parameter in non-robust constraint")

                                        param_id = t.monomial[0][0].id
                                        # coeff = self._float_repr(t.coeff)

                                        if self._prob1A.get(succ):
                                            cons = cons + coeff * self._paramVars[param_id]

                                        # Add the convexified terms
                                        else:
                                            cons = cons + coeff * self._pinit[succ] * (self._paramVars[param_id])
                                            cons = cons + coeff * self._paraminit[param_id] * (
                                                        self._pVars[succ] - self._pinit[succ])

                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if self._prob1A.get(succ):
                                    cons = cons + self._float_repr(constant_value)
                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + self._float_repr(constant_value) * self._pVars[succ]

                        try:
                            # additional step to make the model more stable by removing
                            # small enough constants, they should be accompained by the slacks
                            if abs(cons.getConstant()) < 1e-6:
                                cons = cons - cons.getConstant()
                        except:
                            pass
                        self._remove_set.append(
                            self._encoding.addConstr(self._pVars[state.id] <= cons + self._tau[state.id]))

        end = time.time()


    def _set_objective(self,model):
        """
        Creates the objective
        :param model: stormpy model for uPOMDP
        :return:
        """
        self._objective = 0.0
        # Adding terms to the objective
        for state in model.states:
            self._objective = self._objective + self._tau[state.id]
            self._objective = self._objective - self._pVars[state.id] / self._mu


    def run(self, model, parameters, interval_parameters, properties, prob0E, prob1A, threshold, direction, options, intervals, items,model_check):
        """
        Runs the QCQP procedure by a series of calls to gurobi.

        :param model: The model
        :type model: a stormpy dtmc/mdp
        :param parameters: The policy parameters occuring in the model
        :type parameters: a list of pycarl variables
        :param interval_parameters: The uncertain parameters occuring in the model
        :type parameters: a list of pycarl variables
        :param properties: The properties as an iterable over stormpy.properties
        :param prob0E: The states of having a probability 0 of satisfying the specification for any valuation of
        policies and uncertain parameters
        :type Prob0E: a set
        :param prob1A: The states of having a probability 1 of satisfying the specification for any valuation of
        policies and uncertain parameters
        :type Prob1A: a set
        :param threshold: The threshold
        :type threshold: float
        :param direction: Are we looking for a value below or above
        :type direction: a string, either "above" or "below", we assume above in this code
        :param options: Further options with which the algorithm should run
        :param intervals: List of uncertain intervals for each state and their lower and upper bounds
        :param items: List of uncertain intervals in the model and their lower and upper bounds
        :param model_check: boolean value whether we use robust model checker or not, assumed to be True
        :return:
        """

        #auxillary variables for storing
        self._prob0E = prob0E
        self._prob1A = prob1A
        self._parameters = parameters
        self._interval_parameters= interval_parameters
        self._properties= properties
        self._threshold= threshold
        self._options= options
        self._intervals=intervals
        self._items= items
        self._model_check = model_check
        self._model= model
        self._mu = options.mu
        self._remove_set = []

        assert direction in ["above", "below"]
        if direction == "above":
            raise RuntimeError("Direction == above is currently not supported.")
        if not options.silent:
            print("Number of pmc states: {}".format(model.nr_states))
            print("Number of pmc transitions: {}".format(model.nr_transitions))
            print("Labels: {}".format(model.labeling.get_labels()))
            print(model.model_type)
            print("Number of states: {}".format(model.nr_states))

        results = []

        numstate = model.nr_states
        #Initial solution is uniform over policies
        solution=dict()
        for x in self._parameters:
            solution[x]=stormpy.RationalRF(0.5)
        #gather the values of the interval for robust value checking
        regiondict=dict()
        #print("region check")
        for x in self._interval_parameters:
            for item in items:
                if item.name==x.name:
                    regiondict[x]=(stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))
        #computes the region using stormpy
        region = stormpy.pars.ParameterRegion(regiondict)
        print(region)
        #creates the uncertain MC
        instantiator = stormpy.pars.PartialPDtmcInstantiator(model)
        instantiated_model = instantiator.instantiate(solution)


        ## creates robust verification checker
        env = stormpy.Environment()
        start_check=time.time()
        #robust verification steps
        region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
                                                            allow_model_simplification=False)
        result = region_checker.get_bound_all_states(env, region, maximise=False)
        end_check=time.time()
        self.model_check_timer+=end_check-start_check
        #gather initial solution to the policy parameters for convexification
        self._paraminit = dict([[x.id, float(stormpy.RationalRF(solution[x]))] for x in parameters])
        #print("total model check time:",self.model_check_timer)
        initstate = int(model.initial_states[0])
        #print("model checking ans:")

        #initialize trust region
        trust_region=1.5
        #gather initial solution to the probability variables from robust verification for convexification
        self._pinit = [threshold for _ in range(numstate)]
        for state in model.states:
            self._pinit[int(state.id)]=(result.at(state))
        bestval=self._pinit[initstate]
        #print(bestval)



        # The set of uncertain states
        robuststates1 = []
        for item in intervals:
            stateval = int(item.state)
            # print(stateval,state.id)
            # if stateval==int(state.id):
            succval = int(item.successor)
            robuststates1.append(stateval)
        self._robuststates = set(robuststates1)
        #call the functions to build encoding
        self._create_encoding(model)
        self._set_objective(model)

        # add specification consraint
        self._encoding.addConstr(self._pVars[initstate] + self._tau[initstate] >= self._threshold)
        #main lpp[
        print("-------------\n\nresult, time")
        print("{0}, {1}".format(bestval, self.model_check_timer))
        results.append((bestval, self.model_check_timer))
        for i in range(1,options.maxiter):
            self.iterations = i
            start= time.time()
            # build model constraints
            self._model_constraints(model,i)
            # build trust region constraints
            for x in parameters:
                self._remove_set.append(self._encoding.addConstr(self._paramVars[x.id]<=trust_region* self._paraminit[x.id]))
                self._remove_set.append(self._encoding.addConstr(self._paramVars[x.id]>=self._paraminit[x.id]/trust_region))
            for state in model.states:
                self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] <= trust_region * self._pinit[state.id]))
                self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] >= self._pinit[state.id] / trust_region))
            end = time.time()
            self.encoding_timer += (end - start)

            #solve the problem with Gurobi
            start3 = time.time()
            self._encoding.setObjective(self._objective, GRB.MINIMIZE)
            #print('Solving...')
            self._encoding.optimize()
            t3 = time.time()
            self.solver_timer += (t3 - start3)
            #print("Solver time :" + str(t3 - start3))
            #print("total solver time:",self.solver_timer)
            #print("total encoding time:",self.encoding_timer)

            #print("num iteration",i)

            # Prints the maximum violation
            maxx = 0
            try:
                for state in range(numstate):
                    val = self._tau[state].x
                    if val > maxx:
                        maxx = val
                    #val = tau_neg[state].x
                    #if val > maxx:
                    #    maxx = val

                #if not options.silent:
                    #print("Max vio :", maxx)
                    #print("p =", pVars[initstate].x)

                #gather the parameter values for policy variables
                parameter_values = dict([[id, param_var.x] for id, param_var in self._paramVars.items()])
                #if robust model checking is allowed
                if model_check:
                    #gather the solution for policy parameters
                    solution = dict()
                    for x in self._parameters:
                        solution[x] = stormpy.RationalRF(parameter_values[x.id])
                    #create the region with uncertain parameters
                    regiondict = dict()
                    for x in self._interval_parameters:
                        for item in self._items:
                            if item.name == x.name:
                                regiondict[x] = (stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))
                    # creates the robust model checker with storm
                    region = stormpy.pars.ParameterRegion(regiondict)
                    instantiator = stormpy.pars.PartialPDtmcInstantiator(model)
                    instantiated_model = instantiator.instantiate(solution)
                    env = stormpy.Environment()
                    start_check = time.time()
                    #robust model checking step
                    region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
                                                                        allow_model_simplification=False)

                    #print("region check")
                    result = region_checker.get_bound_all_states(env, region, maximise=False)
                    end_check = time.time()
                    self.model_check_timer += end_check - start_check

                    self._paraminit = dict([[x.id, float(stormpy.RationalRF(solution[x]))] for x in self._parameters])
                    #print("model checking ans:")
                    ansval=result.at(model.initial_states[0])
                    #print(ansval)
                    #if the solution is better than the threshold, return the solution
                    if ansval > threshold:
                        print("Early termination due to positive model checking result at iteration {0}: ".format(str(i)))
                        #print("p[init] = " + str(ansval) )
                        #print("SCP parameter values: ")
                        #print("total model check time:",self.model_check_timer)
                        self.solver_params=solution
                        #print("{0}, {1}".format(ansval, self.model_check_timer+self.solver_timer+self.encoding_timer))
                        results.append((ansval, self.model_check_timer+self.solver_timer+self.encoding_timer))
                        print("{0}, {1}".format(bestval,
                                                self.model_check_timer + self.solver_timer + self.encoding_timer))

                        return results, "policy evaluation passed threshold"

                        #return QcqpResult(self._pVars[initstate].x, parameter_values)
                    #if the solution is better than the best solution so far, update the probability variables and
                    #policy parameters
                    elif ansval>bestval:

                        bestval=ansval
                        results.append((bestval,self.model_check_timer+self.solver_timer+self.encoding_timer))

                        print("{0}, {1}".format(bestval, self.model_check_timer + self.solver_timer + self.encoding_timer))

                        #(self.solver_output.append([bestval,self.model_check_timer+self.solver_timer+self.encoding_timer]))
                        #print("best found values and their computation so far")
                        #print(self.solver_output)
                        #update the probability values for next ieration
                        for state in model.states:
                            self._pinit[int(state.id)] = (result.at(state))
                        #print(pinit)

                        # Updares the policy parameter values for next iteration
                        for param_id, param_var in self._paramVars.items():
                            if not isinstance(param_var, int):
                                self._paraminit[param_id] = param_var.x
                        #increase the size of trust region
                        trust_region=min(10,(trust_region-1)*1.5+1)
                        self.solver_params=solution

                    #shrink the size of trust region
                    else:
                        trust_region=((trust_region-1)/1.5+1)

                    self._encoding.update()

                    #terminate if the trust region is too small
                    if trust_region<1+1e-4:
                        print("Early termination due to small trust region {0}: ".format(str(i)))
                        print("{0}, {1}".format(bestval,
                                                self.model_check_timer + self.solver_timer + self.encoding_timer))

                        return results, "small trust region"
                        #print("p[init] = " + str(ansval))
                        #print("SCP parameter values: ")
                        #for id, param_var in paramVars.items():
                        #    print(str(parameter_names[id]) + "  :  " + str(param_var.x))
                        #return QcqpResult(self._pVars[initstate].x, parameter_values)
                    #print("bestval:",bestval)
            # if gurobi runs into numerical trouble, shrink the trust region
            except AttributeError:
                trust_region = ((trust_region - 1) / 1.5 + 1)
                if trust_region < 1 + 1e-4:
                    print("Early termination due to small trust region {0}: ".format(str(i)))
                    print("{0}, {1}".format(bestval, self.model_check_timer + self.solver_timer + self.encoding_timer))

                    return results, "small trust region"
                    #print("p[init] = " + str(bestval))
                    #print("SCP parameter values: ")
                    #return QcqpResult(bestval, self._paraminit)
                self._encoding.update()
            #print("total model check time:",self.model_check_timer)
            #if the computation time exceeds the time-out, break
            if self.model_check_timer+self.solver_timer+self.encoding_timer>self._options.timeout:
                print("terminating because of timeout")
                print("{0}, {1}".format(bestval, self.model_check_timer + self.solver_timer + self.encoding_timer))
                return results, "timeout of > {}".format(self._options.timeout)
                #print("printing best values and seconds")
                #for item in self.solver_output:
                #    print(item[0], item[1])
                #self.solver_params=solution
                #break
            #remove the constraints that are convexified
            self._encoding.remove(self._remove_set)
            self._remove_set = []

            self._encoding.update()
        print("termination due to max iterations reached: " + str(options.maxiter))
        return results, "max iterations reached"

# function that runs
def test():

    #model path
    path = "pomdp_attempt_prob_36_sat_065_dist_5_obs_diff_orb_len_1.drn"
    #the specification string
    formula_str = "P=? [F \"goal\"]"
    # building options for stormpy
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    # build the model
    pomdp = stormpy.build_parametric_model_from_drn(path, opts)
    # path for the interval values
    interval_path="satellite_robust.intervals"
    # gather the interval parameters
    pomdp_parameters = pomdp.collect_probability_parameters()
    # make the pomdp canonic, which is required for simple transformation.
    pomdp = stormpy.pomdp.make_canonic(pomdp)
    # construct the memory for the FSC
    memory_builder = stormpy.pomdp.PomdpMemoryBuilder()
    #number of memory states in the FSC
    memval=1
    # build the memory
    memory = memory_builder.build(stormpy.pomdp.PomdpMemoryPattern.selective_counter, memval)
    # apply the memory onto the POMDP to get the cartesian product
    pomdp = stormpy.pomdp.unfold_memory(pomdp, memory)
    print("Number of pomdp states before simple:",pomdp.nr_states)
    print("Number of transitions: {}".format(pomdp.nr_transitions))

    # make the POMDP simple.
    pomdp = stormpy.pomdp.make_simple(pomdp)
    print("Number of pomdp states after simple:",pomdp.nr_states)
    print("Number of transitions: {}".format(pomdp.nr_transitions))

    # apply the unknown FSC to obtain a pmc from the POMDP
    pmc = stormpy.pomdp.apply_unknown_fsc(pomdp, stormpy.pomdp.PomdpFscApplicationMode.simple_linear)
    print("Number of pomdp states after simple",pmc.nr_states)
    print("Number of transitions: {}".format(pmc.nr_transitions))
    print("applied pmc")
    path_pmc = "export_" + str(memval) + "_mem"
    #export the resulting model file
    stormpy.export_parametric_to_drn(pmc, path_pmc)
    print("built model")
    fsc_parameters = pmc.collect_probability_parameters() - pomdp.collect_probability_parameters()
    print("number of pomdp parameters:",len(fsc_parameters))

    # gather intervals and items from the model
    intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)

    for item in items:
        print(item,"printing in the main")

    properties = stormpy.parse_properties(formula_str)
    print("Building model from {}".format(path))
    print(pmc.nr_states)
    #compute prob01max states for any policy
    prob0E, prob1A = stormpy.prob01max_states(pmc, properties[0].raw_formula.subformula)
    #compute threshold
    threshold = 0.9999
    direction = "below"  # can be "below" or "above"
    # solver parameters
    options = QcqpOptions(mu=1e4, maxiter=200, graph_epsilon=1e-2, silent=False)

    # run solver
    solver = QcqpSolver_affine_simple_fun()
    result = solver.run(pmc,  fsc_parameters, pomdp_parameters,properties, prob0E, prob1A, threshold, direction, options,intervals,items,True)

    print("number of iterations={}".format(solver.iterations))
    print("solver time={}".format(solver.solver_timer))

    print("solver time={}".format(solver.solver_timer))


def main_prism(path, interval_path, formula_str, threshold, memval=1,timeout=1800, maxiter=200, evaluation_set = []):
    prism_program = stormpy.parse_prism_program(path)

    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    properties = stormpy.parse_properties_for_prism_program(formula_str, prism_program)

    pomdp = stormpy.build_parametric_model(prism_program, properties)
    pomdp_parameters = pomdp.collect_probability_parameters()
    return main(pomdp, interval_path, formula_str, threshold, memval, path, timeout, maxiter, evaluation_set)


def main_drn(path, interval_path, formula_str, threshold, memval=1,timeout=1800, maxiter=200, evaluation_set = []):
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    pomdp = stormpy.build_parametric_model_from_drn(path, opts)
    return main(pomdp, interval_path, formula_str, threshold, memval, path, timeout, maxiter, evaluation_set)


def main(pomdp, interval_path, formula_str, threshold, memval=1, path="",timeout=1800, maxiter=200, evaluation_set = []):


    model_info = dict()
    model_info["model name"] = path
    model_info["interval file"] = interval_path
    model_info["objective"] = formula_str
    model_info["timeout"] = timeout
    model_info["max iterations"] = maxiter

    t0 = time.time()

    pomdp_parameters = pomdp.collect_probability_parameters()
    #stormpy.export_to_drn(pomdp, "pomdp_ex")
    pomdp = stormpy.pomdp.make_canonic(pomdp)

    model_info["nr_model_states"] = pomdp.nr_states
    model_info["nr_model_transition"] = pomdp.nr_transitions

    memory_builder = stormpy.pomdp.PomdpMemoryBuilder()
    memory = memory_builder.build(stormpy.pomdp.PomdpMemoryPattern.selective_counter, memval)
    pomdp = stormpy.pomdp.unfold_memory(pomdp, memory)
    pomdp = stormpy.pomdp.make_simple(pomdp)

    model_info["nr_product_states"] = pomdp.nr_states
    model_info["nr_product_transitions"] = pomdp.nr_transitions
    model_info["mem_size"] = memval

    pmc = stormpy.pomdp.apply_unknown_fsc(pomdp, stormpy.pomdp.PomdpFscApplicationMode.simple_linear)


    fsc_parameters = pmc.collect_probability_parameters() - pomdp.collect_probability_parameters()

    intervals, polyhedrons, items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)


    properties = stormpy.parse_properties(formula_str)
    prob0E, prob1A = stormpy.prob01max_states(pmc, properties[0].raw_formula.subformula)
    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=1e4, maxiter=maxiter, graph_epsilon=1e-2, silent=False,timeout=timeout)

    print("Model info:")
    print(model_info)
    print("\n")

    t1 = time.time()
    solver = QcqpSolver_affine_simple_fun()
    solver_results, solver_exit = solver.run(pmc,  fsc_parameters, pomdp_parameters,properties, prob0E, prob1A, threshold, direction, options,intervals,items,True)

    model_info["total iterations"] = solver.iterations
    model_info["total solver time"] = solver.solver_timer
    model_info["exit code"] = solver_exit




    interval_set = [interval_path] + evaluation_set
    eval_results = dict()
    for eval_set_path in interval_set:
        intervals, polyhedrons, items = interval_parser.parse_model_interval(pmc, pomdp_parameters, eval_set_path)

        regiondict = dict()
        for x in pomdp_parameters:
            for item in items:
                if item.name == x.name:
                    regiondict[x] = (stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))
        region = stormpy.pars.ParameterRegion(regiondict)
        instantiator = stormpy.pars.PartialPDtmcInstantiator(pmc)
        instantiated_model = instantiator.instantiate(solver.solver_params)
        env = stormpy.Environment()
        env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
        env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
        env.solver_environment.native_solver_environment.precision = stormpy.Rational('0.01')
        region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
                                                            allow_model_simplification=False)
        result = region_checker.get_bound_all_states(env, region, maximise=False)
        eval_results[eval_set_path] = result.at(pmc.initial_states[0])
        #ansval = result.at(pmc.initial_states[0])
        if eval_set_path == interval_path:
            t1_end = time.time()
            solver_results.append((eval_results[interval_path], t1_end - t1))

    tend = time.time()

    print("Final result: ", eval_results[interval_path])
    print("Total solver time: ", str(tend - t1))
    print("Total computation time: ", str(tend - t0))

    model_info["total solver time"] = str(tend - t1)
    model_info["total computation time"] = str(tend - t0)
    model_info["final result"] = eval_results[interval_path]

    return model_info, solver_results, eval_results






if __name__ == '__main__':
    path, interval_path, formula_str, threshold = input_files()
    main(path, interval_path, formula_str, threshold)
