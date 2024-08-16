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


class QcqpOptions():
    def __init__(self, mu, maxiter, graph_epsilon, silent,timeout):
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
        Returns the float representation for a constant value
        :model: stormpy model of uPOMDP
        :return:
        """
        numstate = model.nr_states

        self._encoding = Model("qcp")
        self._encoding.setParam('OutputFlag', not self._options.silent)

        # Initializing some arrays for state, parameter and tau variables, and their values at previous iterations
        self._pVars = [self._encoding.addVar(lb=0) for _ in range(numstate)]
        #maxact_list = [0 for _ in range(numstate)]
        self._tau = [self._encoding.addVar(lb=0) for _ in range(numstate)]
        # tau_neg = [m.addVar(lb=0) for _ in range(numstate)]
        # Initializing gurobi variables for parameters,lb=lowerbound, ub=upperbound
        self._paramVars = dict([[x.id, self._encoding.addVar(lb=self._options.graph_epsilon, ub=1 - self._options.graph_epsilon)] for x in self._parameters])
        #dicts for dual variables of the semi-infinite constraints
        self._dual_upper = dict()
        self._dual_lower = dict()
        for item in self._intervals:
            stateval = int(item.state)
            # print(stateval,state.id)
            # if stateval==int(state.id):
            succval = int(item.successor)
            self._dual_upper[(stateval, succval)] = self._encoding.addVar(lb=0)
            self._dual_lower[(stateval, succval)] = self._encoding.addVar(lb=0)
        self._encoding.update()
        #gurobi params
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
        Returns the float representation for a constant value
        :model: stormpy model of uPOMDP
        :i: current iteration, it is important if it is 1 or not
        :return:
        """

        for state in model.states:
            #print the number of states in the loop
            if int(state.id) % 10000 == 0:
                print("Encoding through states",state.id)
            start2 = time.time()
            #for initial iteration, build everything, including uncertain constraints
            if i==1:
                #if not in robust states
                if int(state.id) not in self._robuststates:
                    start2 = time.time()

                    # go over actions, usually 1 for uPOMDPs
                    for action in state.actions:
                        reward_at_state=self._reward_name[state.id]
                        #print(reward_at_state)
                        if reward_at_state.is_constant():
                            cons = self._float_repr(reward_at_state.constant_part())
                        else:
                            cons=0
                            num = reward_at_state.numerator.polynomial()
                            den = reward_at_state.denominator.constant_part()

                            for t in num:
                                coeff = self._float_repr(t.coeff / den)
                                if t.is_constant():
                                    cons=cons+coeff
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    # print(t.monomial[0][0],state)
                                    if t.monomial[0][0] in self._interval_parameters:
                                        raise RuntimeError("Interval parameter in non-robust constraint")

                                    param_id = t.monomial[0][0].id
                                    cons = cons + coeff  * (self._paramVars[param_id])

                            #print(reward_at_state,state.id)
                        #go over transitions
                        for transition in action.transitions:

                            succ = int(transition.column)
                            # Value of transition
                            transition_value = transition.value()
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
                                        if self._rew0E.get(succ):
                                            pass
                                        # Add nothing successor state is prob0
                                        # elif prob0E.get(succ):
                                        #    pass
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            cons = cons + coeff * self._pVars[succ]

                                    # If the transition term is not a constant
                                    else:
                                        if t.tdeg > 1:
                                            raise RuntimeError("We expect the term to be a single variable")
                                        # print(t.monomial[0][0],state)
                                        if t.monomial[0][0] in self._interval_parameters:
                                            raise RuntimeError("Interval parameter in non-robust constraint")

                                        param_id = t.monomial[0][0].id
                                        # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                        if  self._rew0E.get(succ):
                                            pass
                                            #cons = cons + coeff * self._paramVars[param_id]
                                        # Add nothing successor state is prob0
                                        # elif prob0E.get(succ):
                                        #    pass
                                        # Add the convexified terms, this is now affine
                                        else:
                                            cons = cons + coeff * self._pinit[succ] * (self._paramVars[param_id])
                                            cons = cons + coeff * self._paraminit[param_id] * (self._pVars[succ] - self._pinit[succ])

                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if self._rew0E.get(succ):
                                    pass
                                # If successor state is prob0, do nothing
                                # elif prob0E.get(succ):
                                #    pass
                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + self._float_repr(constant_value) * self._pVars[succ]

                        try:
                            # additional step to make the model more stable by removing
                            # small enough constants, they should be accompained by the slacks
                            if abs(cons.getConstant()) < 1e-4:
                                cons = cons - cons.getConstant()
                        except:
                            pass
                        #add these constraints to "remove_set", which we remove after each solver iteration
                        self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] >= cons - self._tau[state.id]))
                    end2 = time.time()
                    self.robust_encoding_timer += (end2 - start2)
                #robust constraints
                else:
                    start2 = time.time()

                    for action in state.actions:
                        reward_at_state=self._reward_name[state.id]
                        #print(reward_at_state)
                        if reward_at_state.is_constant():
                            cons = self._float_repr(reward_at_state.constant_part())
                        else:
                            cons=0
                            num = reward_at_state.numerator.polynomial()
                            den = reward_at_state.denominator.constant_part()

                            for t in num:
                                coeff = self._float_repr(t.coeff / den)
                                if t.is_constant():
                                    cons=cons+coeff
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    # print(t.monomial[0][0],state)
                                    if t.monomial[0][0] in self._interval_parameters:
                                        raise RuntimeError("Interval parameter in non-robust constraint")

                                    param_id = t.monomial[0][0].id
                                    cons = cons + coeff  * (self._paramVars[param_id])

                            #print(reward_at_state,state.id)
                        for transition in action.transitions:

                            succ = int(transition.column)
                            # Value of transition
                            transition_value = transition.value()
                            # TODO check that the denominator is indeed constant?
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
                                        if self._rew0E.get(succ):
                                            pass
                                        # Add nothing successor state is prob0
                                        # elif prob0E.get(succ):
                                        #    pass
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            cons = cons + coeff * self._pVars[succ]

                                    # If the transition term is not a constant
                                    else:
                                        if t.tdeg == 1:
                                            # if the parameter is in interval

                                            if t.monomial[0][0] in self._interval_parameters:
                                                interval_name = t.monomial[0][0].name
                                                interval_value = -1
                                                #get the name of param
                                                for item in self._items:
                                                    if item.name == interval_name:
                                                        interval_value = 0
                                                        break
                                                if interval_value == -1:
                                                    # print("interval name: " + str(interval_name))
                                                    raise RuntimeError("Interval not part of the instantiation")
                                                stateval = int(state.id)
                                                succval = int(succ)
                                                # add robust constraint, and terms to the initial constraint due to dualization
                                                self._encoding.addConstr(self._dual_upper[(stateval, succval)] - self._dual_lower[
                                                    (stateval, succval)] == coeff * self._pVars[succ])
                                                cons = cons +self._dual_upper[(stateval, succval)] * item.upperbound - \
                                                       self._dual_lower[(stateval, succval)] * item.lowerbound
                                                # cons= cons + coeff*item.upperbound*pVars[succ]
                                            # we should not reach this part for simple pomdps, adding here for transitions with
                                            # fsc and uncertain parameters
                                            elif t.monomial[0][0] in self._parameters:
                                                raise RuntimeError("Assumed the POMDP is simple")


                                        elif t.tdeg > 1:
                                            raise RuntimeError("the supplied POMDP is not simple")
                            # If the value of transition is constant
                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if rew0E.get(succ):
                                    pass
                                # If successor state is prob0, do nothing
                                # elif prob0E.get(succ):
                                #    pass
                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + self._float_repr(constant_value) * self._pVars[succ]
                        try:
                            # additional step to make the model more stable by removing
                            # small enough constants, they should be accompained by the slacks
                            if abs(cons.getConstant()) < 1e-4:
                                cons = cons - cons.getConstant()
                        except:
                            pass
                        #do not add to remove
                        self._encoding.addConstr(self._pVars[state.id] >= cons - self._tau[state.id])
                    end2 = time.time()
                    self.robust_encoding_timer += (end2 - start2)
            else:
                #just go over the nonrobust states for sucessive iterations
                if int(state.id) not in self._robuststates:
                    start2 = time.time()

                    # non robust constraint:
                    for action in state.actions:
                        reward_at_state=self._reward_name[state.id]
                        #print(reward_at_state)
                        if reward_at_state.is_constant():
                            cons = self._float_repr(reward_at_state.constant_part())
                        else:
                            cons=0
                            num = reward_at_state.numerator.polynomial()
                            den = reward_at_state.denominator.constant_part()

                            for t in num:
                                coeff = self._float_repr(t.coeff / den)
                                if t.is_constant():
                                    cons=cons+coeff
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    # print(t.monomial[0][0],state)
                                    if t.monomial[0][0] in self._interval_parameters:
                                        raise RuntimeError("Interval parameter in non-robust constraint")

                                    param_id = t.monomial[0][0].id
                                    cons = cons + coeff  * (self._paramVars[param_id])

                            #print(reward_at_state,state.id)
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
                                        if self._rew0E.get(succ):
                                            pass
                                        # Add nothing successor state is prob0
                                        # elif prob0E.get(succ):
                                        #    pass
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            cons = cons + coeff * self._pVars[succ]

                                    # If the transition term is not a constant
                                    else:
                                        # print(t)
                                        if t.tdeg > 1:
                                            raise RuntimeError("We expect the term to be a single variable")
                                        # print(t.monomial[0][0],state)
                                        if t.monomial[0][0] in self._interval_parameters:
                                            raise RuntimeError("Interval parameter in non-robust constraint")

                                        param_id = t.monomial[0][0].id
                                        # coeff = self._float_repr(t.coeff)

                                        # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                        if self._rew0E.get(succ):
                                            cons = cons + coeff * self._paramVars[param_id]
                                        # Add nothing successor state is prob0
                                        # elif prob0E.get(succ):
                                        #    pass
                                        # Add the convexified terms
                                        else:
                                            cons = cons + coeff * self._pinit[succ] * (self._paramVars[param_id])
                                            cons = cons + coeff * self._paraminit[param_id] * (
                                                        self._pVars[succ] - self._pinit[succ])

                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if self._rew0E.get(succ):
                                    pass
                                # If successor state is prob0, do nothing
                                # elif prob0E.get(succ):
                                #    pass
                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + self._float_repr(constant_value) * self._pVars[succ]

                        try:
                            # additional step to make the model more stable by removing
                            # small enough constants, they should be accompained by the slacks
                            if abs(cons.getConstant()) < 1e-4:
                                cons = cons - cons.getConstant()
                        except:
                            pass
                        self._remove_set.append(
                            self._encoding.addConstr(self._pVars[state.id] >= cons - self._tau[state.id]))
        # for state in model.states:
        #     self._remove_set.append(
        #         self._encoding.addConstr(self._pVars[state.id] <= trust_region * self._pinit[state.id]))
        #     self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] >= self._pinit[state.id] / trust_region))

        end = time.time()

    def _set_objective(self,model):
        """
        Sets the objective of the optimization problem
        :model: stormpy model of uPOMDP
        :return:
        """

        self._objective = 0.0
        # Adding terms to the objective
        for state in model.states:
            self._objective = self._objective + self._tau[state.id]
            # objective = objective + tau_neg[state.id]
            self._objective = self._objective + self._pVars[state.id] / self._mu


    def run(self, model, parameters, interval_parameters, properties, rew0E, reward_model_name, threshold, direction, options, intervals,items,model_check):
        """
        Runs the SCP procedure by a series of calls to gurobi.

        :param model: The model
        :type model: a stormpy dtmc/mdp
        :param parameters: The FSC parameters occuring in the model
        :type parameters: a list of pycarl variables
        :param interval_parameters: The uncertain parameters occuring in the model
        :type interval_parameters: a list of pycarl variables
        :param properties: The properties as an iterable over stormpy.properties
        :param rew0E: set of states with rew0
        :param reward_model_name: name of the reward model
        :param threshold: The threshold
        :type threshold: float
        :param direction: Are we looking for a value below or above
        :type direction: a string, either "above" or "below"
        :param options: Further options with which the algorithm should run
        :param items: Set of intervals with lower and upper value
        :param model_check: boolean value if we model check at each iteration
        :return:
        """

        #storing items
        self._rew0E = rew0E
        #self._prob1A = prob1A
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

        reward_name = model.reward_models[reward_model_name]
        self._reward_name=reward_name.state_rewards
        #self._reward_dict = dict()

        numstate = model.nr_states
        # print(numstate)
        # Initializing all params for 0.5, and perform a model checking step
        solution=dict()
        for x in self._parameters:
            solution[x]=stormpy.RationalRF(0.5)
        regiondict=dict()
        print("region check")
        # Go over interval params to create regions
        for x in self._interval_parameters:
            for item in items:
                if item.name==x.name:
                    regiondict[x]=(stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))

        #Instantiate the robust model checking
        region = stormpy.pars.ParameterRegion(regiondict)
        print(region)
        instantiator = stormpy.pars.PartialPDtmcInstantiator(model)
        instantiated_model = instantiator.instantiate(solution)

        ## create epsilon-region
        env = stormpy.Environment()
        env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
        env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration

        start_check=time.time()
        region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
                                                            allow_model_simplification=False)
        #get results for each state
        result = region_checker.get_bound_all_states(env, region, maximise=True)
        end_check=time.time()

        self.model_check_timer+=end_check-start_check
        #initialize param values
        self._paraminit = dict([[x.id, float(stormpy.RationalRF(solution[x]))] for x in parameters])
        print("total model check time:",self.model_check_timer)
        #print(paraminit)
        initstate = int(model.initial_states[0])
        #
        #initialsize size of trust region
        trust_region=1.5
        #initialize probability values
        self._pinit = [0 for _ in range(numstate)]
        for state in model.states:
            self._pinit[int(state.id)]=(result.at(state))
        #print(pinit)
        #set best found value
        bestval=self._pinit[initstate]
        print("model checking ans:",self._pinit[initstate])

        results.append((bestval, self.model_check_timer))
        #print(bestval)

        # The penalty parameter for constraint violation

        # Compute the set of states with robust transition
        robuststates_set = []
        # Updates the model for gurobi
        # for state in model.states:
        # maxact=0
        for item in intervals:
            stateval = int(item.state)
            # print(stateval,state.id)
            # if stateval==int(state.id):
            succval = int(item.successor)
            robuststates_set.append(stateval)
        self._robuststates = set(robuststates_set)
        #build encoding and objective
        self._create_encoding(model)
        self._set_objective(model)

        #threshold constraint
        self._encoding.addConstr(self._pVars[initstate] <= self._threshold+self._tau[state.id])

        #trust region constraints for  param variables
        for i in range(1,options.maxiter):
            self.iterations = i
            start= time.time()
            self._model_constraints(model,i)

            for x in parameters:
                self._remove_set.append(self._encoding.addConstr(self._paramVars[x.id]<=trust_region* self._paraminit[x.id]))
                self._remove_set.append(self._encoding.addConstr(self._paramVars[x.id]>=self._paraminit[x.id]/trust_region))
            for state in model.states:
                if abs(self._pinit[state.id])>1e-4:
                    self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] <= trust_region * self._pinit[state.id]))
                    self._remove_set.append(self._encoding.addConstr(self._pVars[state.id] >= self._pinit[state.id] / trust_region))

            end = time.time()
            self.encoding_timer += (end - start)

            start3 = time.time()
            self._encoding.setObjective(self._objective, GRB.MINIMIZE)
            print('Solving...')
            self._encoding.optimize()
            t3 = time.time()
            self.solver_timer += (t3 - start3)
            print("Solver time :" + str(t3 - start3))
            print("total solver time:",self.solver_timer)
            print("total encoding time:",self.encoding_timer)
            print("total robust encoding time:",self.robust_encoding_timer)

            print("num iteration",i)
            print("trust_region",trust_region)


            # Prints the maximum violation
            maxx = 0
            try:
                for state in range(numstate):
                    val = self._tau[state].x
                    if val > maxx:
                        maxx = val

                if not options.silent:
                    print("Max vio :", maxx)
                    #print("p =", pVars[initstate].x)


                parameter_values = dict([[id, param_var.x] for id, param_var in self._paramVars.items()])

                #model checking step
                if model_check:
                    #print(solution)
                    option = "max" # min or max for value iterator
                    solution = dict()
                    for x in self._parameters:
                        solution[x] = stormpy.RationalRF(parameter_values[x.id])
                    #print(solution)
                    regiondict = dict()
                    for x in self._interval_parameters:
                        for item in self._items:
                            if item.name == x.name:
                                regiondict[x] = (stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))
                    region = stormpy.pars.ParameterRegion(regiondict)
                    instantiator = stormpy.pars.PartialPDtmcInstantiator(model)
                    instantiated_model = instantiator.instantiate(solution)
                    #assert interval_parameters == instantiated_model.collect_probability_parameters()

                    env = stormpy.Environment()
                    #env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
                    #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration
                    start_check = time.time()

                    region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
                                                                        allow_model_simplification=False)

                    print("region check")
                    result = region_checker.get_bound_all_states(env, region, maximise=True)
                    end_check = time.time()
                    self.model_check_timer += end_check - start_check

                    self._paraminit = dict([[x.id, float(stormpy.RationalRF(solution[x]))] for x in self._parameters])
                    print("model checking ans:")
                    ansval=result.at(model.initial_states[0])
                    print(ansval)
                    #terminate if the result is better than threshold
                    if ansval < threshold:

                        # TODO adjust result
                        print("Early termination due to positive model checking result at iteration {0}: ".format(str(i)))
                        #print("p[init] = " + str(ansval) )
                        #print("SCP parameter values: ")
                        #for id, param_var in paramVars.items():
                        #    print(str(parameter_names[id]) + "  :  " + str(param_var.x))
                        #print("total model check time:",self.model_check_timer)
                        self.solver_params=solution
                        results.append((ansval, self.model_check_timer+self.solver_timer+self.encoding_timer))
                        return results
                        #return QcqpResult(self._pVars[initstate].x, parameter_values)
                    #enlarge trust region and update probability and parameter values
                    #if the best found solution is better

                    elif ansval<bestval:
                        bestval=ansval
                        #self.solver_output[]=self.model_check_timer+self.solver_timer+self.encoding_timer
                        results.append((bestval, self.model_check_timer+self.solver_timer+self.encoding_timer))
                        #(self.solver_output.append([bestval,self.model_check_timer+self.solver_timer+self.encoding_timer]))
                        #print(self.solver_output)


                        for state in model.states:
                            self._pinit[int(state.id)] = (result.at(state))
                        #print(pinit)

                        # Updares the parameter values for next iteration
                        for param_id, param_var in self._paramVars.items():
                            if not isinstance(param_var, int):
                                if abs(param_var.x) > options.graph_epsilon:
                                    #  print pVar
                                    self._paraminit[param_id] = param_var.x
                                else:
                                    self._paraminit[param_id] = self._options.graph_epsilon
                        trust_region=min(10,(trust_region-1)*1.5+1)
                        self.solver_params=solution

                    #contract trust region if solution is not improved
                    else:
                        trust_region=((trust_region-1)/1.5+1)
                    self._encoding.update()
                    #terminate if trust region is small

                    if trust_region<1+1e-4:
                        print("Early termination due to small trust region {0}: ".format(str(i)))
                        return results
                        #print("p[init] = " + str(ansval))
                        #print("SCP parameter values: ")
                        #for id, param_var in paramVars.items():
                        #    print(str(parameter_names[id]) + "  :  " + str(param_var.x))
                        #return QcqpResult(self._pVars[initstate].x, parameter_values)
                    print("bestval:",bestval)
            #safeguard for gurobi failure, contract trust region if this is the case

            except AttributeError:
                trust_region = ((trust_region -1 ) / 1.5 +1)
                if trust_region <  1+1e-4:
                    print("Early termination due to small trust region {0}: ".format(str(i)))
                    #print("p[init] = " + str(bestval))
                    #print("SCP parameter values: ")
                    return results
                    #return QcqpResult(bestval, self._paraminit)
                self._encoding.update()
            print("total model check time:",self.model_check_timer)
            if self.model_check_timer+self.solver_timer+self.encoding_timer>self._options.timeout:
                print("Termination due to timeout")
                return results
                #for item in self.solver_output:
                #    print(item[0], item[1])
                #self.solver_params=solution
                break
            #remove constraints that are required to be updated
            self._encoding.remove(self._remove_set)
            self._remove_set = []

            self._encoding.update()
        return results

# testing function

def main_drn(path, interval_path, formula_str, threshold, memval=1):
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    pomdp = stormpy.build_parametric_model_from_drn(path, opts)
    return main(pomdp, interval_path, formula_str, threshold, memval, path)


def main_prism(path, interval_path, formula_str, threshold, memval=1):
    prism_program = stormpy.parse_prism_program(path)
    # formula_str = "P=? [!\"bad\" U \"goal\"]"
    # formula_str = "P=? [F \"goal\"]"
    # interval_path="collision_partial_obs_2d_upd_hobs_20_small.intervals"
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    properties = stormpy.parse_properties_for_prism_program(formula_str, prism_program)
    # construct the pPOMDP
    import inspect
    print(inspect.getfullargspec(stormpy.build_parametric_model))
    pomdp = stormpy.build_parametric_model(prism_program, properties)

    return main(pomdp, interval_path, formula_str, threshold, memval, path)


def main(pomdp, interval_path, formula_str, threshold, memval=1, path=""):
    model_info = dict()
    model_info["model name"] = path
    model_info["objective"] = formula_str

    t0 = time.time()

    pomdp_parameters = pomdp.collect_probability_parameters()
    stormpy.export_to_drn(pomdp, "pomdp_ex")
    pomdp = stormpy.pomdp.make_canonic(pomdp)
    memory_builder = stormpy.pomdp.PomdpMemoryBuilder()
    memory = memory_builder.build(stormpy.pomdp.PomdpMemoryPattern.selective_counter, memval)
    pomdp = stormpy.pomdp.unfold_memory(pomdp, memory)
    pomdp = stormpy.pomdp.make_simple(pomdp)

    model_info["nr_model_states"] = pomdp.nr_states
    model_info["nr_model_transition"] = pomdp.nr_transitions

    pmc = stormpy.pomdp.apply_unknown_fsc(pomdp, stormpy.pomdp.PomdpFscApplicationMode.simple_linear)

    model_info["nr_product_states"] = pmc.nr_states
    model_info["nr_product_transitions"] = pmc.nr_transitions
    model_info["mem_size"] = memval

    path_pmc = "export_" + str(memval) + "_mem"

    stormpy.export_to_drn(pmc, path_pmc)

    fsc_parameters = pmc.collect_probability_parameters() - pomdp.collect_probability_parameters()

    intervals, polyhedrons, items = interval_parser.parse_model_interval(pmc, pomdp_parameters, interval_path)
    # for item in intervals:
    # print(item,"printing in the main")
    # for item in items:
    #    print(item,"printing in the main")

    # for p in polyhedrons:
    # print(p)
    #    p.compute_vertices()



    # run solver
    properties = stormpy.parse_properties(formula_str)
    rew0E = find_rew0_states(properties[0], pmc)
    reward_name = list(pmc.reward_models.keys())[0]

    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=1e4, maxiter=10000, graph_epsilon=1e-2, silent=False,timeout=35)

    solver = QcqpSolver_affine_simple_fun()
    solver_results = solver.run(pmc, fsc_parameters, pomdp_parameters, properties, rew0E, reward_name, threshold, direction,
                        options, intervals, items, True)


    model_info["total iterations"] = solver.iterations
    model_info["total solver time"] = solver.solver_timer

    # compute the policy against the robust interval
    # interval_path="collision_partial_obs_2d_upd_hobs_20_small.intervals"
    # intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)
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
    ansval = result.at(pmc.initial_states[0])
    tend = time.time()

    print("Final result: ", ansval)
    print("Total time: ", str(tend - t0))

    model_info["total computation time"] = str(tend - t0)
    model_info["final result"] = ansval

    return model_info, solver_results
















    #construct the pomdp from pomdp-drn
    #path = "pomdp_attempt_prob_rew_36_sat_065_dist_5_obs_diff_orb_len_1.drn"

    #formula_str = "R=? [F \"goal\"]"
    # adapt into lower and upper model


    #interval_path="satellite_rew_robust.intervals"
    # construct the pPOMDP
    #pomdp = stormpy.build_parametric_model(prism_program, properties)

    # pomdp_parameters = pomdp.collect_probability_parameters()
    # make its representation canonic.
    # pomdp = stormpy.pomdp.make_canonic(pomdp)
    #stormpy.export_parametric_to_drn(pomdp, "export_pomdp.drn")

    # construct the memory for the FSC
    # in this case, a selective counter with two states
    # memory_builder = stormpy.pomdp.PomdpMemoryBuilder()
    #number of memory states
    #memval=1
    #memory = memory_builder.build(stormpy.pomdp.PomdpMemoryPattern.selective_counter, memval)
    # apply the memory onto the POMDP to get the cartesian product
    #pomdp = stormpy.pomdp.unfold_memory(pomdp, memory)
    #print("Number of pomdp states before simple:",pomdp.nr_states)
    #print("Number of transitions: {}".format(pomdp.nr_transitions))

    # make the POMDP simple. This step is optional but often beneficial
    #pomdp = stormpy.pomdp.make_simple(pomdp)
    #print("Number of pomdp states after simple:",pomdp.nr_states)
    #print("Number of transitions: {}".format(pomdp.nr_transitions))

    # apply the unknown FSC to obtain a pmc from the POMDP
    #pmc = stormpy.pomdp.apply_unknown_fsc(pomdp, stormpy.pomdp.PomdpFscApplicationMode.simple_linear)
    #print("Number of pomdp states after simple",pmc.nr_states)
    #print("Number of transitions: {}".format(pmc.nr_transitions))
    #print("applied pmc")
    #path_pmc = "export_" + str(memval) + "_mem"
    #stormpy.export_to_drn(pmc, path_pmc)
    #print("built model")
    #print(model.initial_states)
    #fsc_parameters = pmc.collect_probability_parameters() - pomdp.collect_probability_parameters()
    #print("number of pomdp parameters:",len(fsc_parameters))
    #print(fsc_parameters)
    #print(pomdp_parameters)
    # gather intervals and items
    #intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)

    #for item in items:
    #    print(item,"printing in the main")

    #properties = stormpy.parse_properties(formula_str)
    #print("Building model")
    #parameters = model.collect_probability_parameters()
    #print(pmc.nr_states)
    #compute prob01max states
    #rew0E = find_rew0_states(properties[0], pmc)
    #reward_name = list(pmc.reward_models.keys())[0]
    #threshold = 1
    #direction = "below"  # can be "below" or "above"
    #options = QcqpOptions(mu=1e4, maxiter=1000, graph_epsilon=2e-2, silent=False)

    # run solver
    #solver = QcqpSolver_affine_simple_fun()
    #result = solver.run(pmc,  fsc_parameters, pomdp_parameters,properties, rew0E, reward_name,threshold, direction, options,intervals,items,True)

    #print("number of iterations={}".format(solver.iterations))
    #print("solver time={}".format(solver.solver_timer))

    #print("solver time={}".format(solver.solver_timer))
    #compute the policy against the robust interval

    #intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)
    #regiondict = dict()
    #for x in pomdp_parameters:
    #    for item in items:
    #        if item.name == x.name:
    #            regiondict[x] = (stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))
    #region = stormpy.pars.ParameterRegion(regiondict)
    #instantiator = stormpy.pars.PartialPDtmcInstantiator(pmc)
    #instantiated_model = instantiator.instantiate(solver.solver_params)
                    #assert interval_parameters == instantiated_model.collect_probability_parameters()

    #env = stormpy.Environment()
    #env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
    #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration
                    #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.optimistic_value_iteration

    #start_check = time.time()

    #region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
    #                                                                    allow_model_simplification=False)

    #print("region check")
    #result = region_checker.get_bound_all_states(env, region, maximise=False)
    #end_check = time.time()
    #print("model checking ans:")
    #ansval=result.at(pmc.initial_states[0])
    #print(ansval)

def test_drn():
    path = "parametric_maze.pm"
    prism_program = stormpy.parse_prism_program(path)

    #formula_str = "P=? [!\"bad\" U \"goal\"]"
    formula_str = "P=? [F \"goal\"]"
    interval_path="parametric_maze.intervals"

    properties = stormpy.parse_properties_for_prism_program(formula_str, prism_program)
    # construct the pPOMDP
    pomdp = stormpy.build_parametric_model(prism_program, properties)


    path = "export_5_mem_pomdp_attempt_prob_36_sat_065_dist_5_obs_diff_orb_len_1.drn"
    print(path)
    formula_str = "P=? [F \"goal\"]"
    # adapt into lower and upper model
    #opts = stormpy.DirectEncodingParserOptions()
    #opts.build_choice_labels = True
    pmc = stormpy.build_parametric_model_from_drn(path)
    print("number of pmc states ",pmc.nr_states)
    print("Number of transitions: {}".format(pmc.nr_transitions))
    path = "pomdp_attempt_prob_36_sat_065_dist_5_obs_diff_orb_len_1.drn"
    #pomdp = stormpy.build_parametric_model_from_drn(path)
    #print("number of pomdp states",pomdp.nr_states)
    #print("Number of transitions: {}".format(pomdp.nr_transitions))
    #
    # properties = stormpy.parse_properties(formula_str)

    #formula_str = "P=? [!\"bad\" U \"goal\"]"

    interval_path="satellite.intervals"
    #properties = stormpy.parse_properties_for_prism_program(formula_str, prism_program)
    # construct the pPOMDP
    #pomdp = stormpy.build_parametric_model(prism_program, properties)

    #pomdp_parameters = pmc.collect_probability_parameters()
    # make its representation canonic.
    #pomdp = stormpy.pomdp.make_canonic(pomdp)
    #stormpy.export_parametric_to_drn(pomdp, "export_pomdp.drn")

    # construct the memory for the FSC
    # in this case, a selective counter with two states

    #fsc_parameters = pmc.collect_probability_parameters() - pomdp.collect_probability_parameters()
    parameters=pmc.collect_probability_parameters()
    pomdp_parameters=set()
    fsc_parameters=set()
    for item in parameters:
        if not item.name[0]=="p":

            pomdp_parameters.add(item)
        else:
            fsc_parameters.add(item)
            #print(item)
    #print(fsc_parameters)
    #print(pomdp_parameters)
    #print(type(fsc_parameters))
    #print(dir(fsc_parameters))

    #print(pomdp.collect_probability_parameters())
    print("number of pomdp parameters:",len(fsc_parameters))
    #intervals2, polyhedrons = interval_parser.parse_input(interval_path)
    intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)
    #for item in intervals:
        #print(item,"printing in the main")
    for item in items:
        print(item,"printing in the main")
    for p in polyhedrons:
        #print(p)
        p.compute_vertices()



    properties = stormpy.parse_properties(formula_str)
    print("Building model from {}".format(path))
    #parameters = model.collect_probability_parameters()
    print(pmc.nr_states)

    prob0E, prob1A = stormpy.prob01max_states(pmc, properties[0].raw_formula.subformula)

    threshold = 0.1
    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=1e4, maxiter=10000, graph_epsilon=1e-2, silent=False)

    # result = solver.run(reward_model_name ,model_rew, parameters_rew, rew0, rew_threshold, direction, options)
    solver = QcqpSolver_affine_simple_fun()
    result = solver.run(pmc, fsc_parameters, pomdp_parameters,properties, prob0E, prob1A, threshold, direction, options,intervals,items,True)

    print("number of iterations={}".format(solver.iterations))
    print("solver time={}".format(solver.solver_timer))

    print("solver time={}".format(solver.solver_timer))



def example_getting_started_06():
    path = "/home/lmc3696/Desktop/research/codes/qcqp_codes/pomdp/pomdp_casestudies/4x4grid_sl_avoid_standard_m2.drn"
    formula_str = 'P=? [F "goal" ]'

    # adapt into lower and upper model
    model = stormpy.build_parametric_model_from_drn(path)

    properties = stormpy.parse_properties(formula_str)

    #  formulas = stormpy.parse_properties_for_prism_program(formula_str, path)
  #  properties = stormpy.parse_properties_for_prism_program(formula_str, path)
   # model = stormpy.build_parametric_model(path, properties)
    print("Building model from {}".format(path))


    threshold = 0.975
    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=100, maxiter=10000, graph_epsilon=1e-6, silent=False)

    parameters = model.collect_probability_parameters()
    prob0E, prob1A = stormpy.prob01max_states(model, properties[0].raw_formula.subformula)
    print(prob1A)
    #

    #result = solver.run(reward_model_name ,model_rew, parameters_rew, rew0, rew_threshold, direction, options)
    solver = QcqpSolver_affine()
    result = solver.run(model, parameters, properties, prob0E, prob1A, threshold, direction, options, polyhedrons,False)
    parameter_assignments = dict([[x, result.parameter_values[x.id]] for x in parameters])
    #print(parameter_assignments)
    instantiator = stormpy.pars.PDtmcInstantiator(model)

    # # Check distance to result by storm (notice that also the point that storm checks slightly differs)
    rational_parameter_assignments = dict([[x, stormpy.RationalRF(val)] for x, val in parameter_assignments.items()])
    #print(rational_parameter_assignments)
    instantiated_model = instantiator.instantiate(rational_parameter_assignments)
    mc_res = stormpy.model_checking(instantiated_model, properties[0]).at(model.initial_states[0])

    print("Qcqp: {}: Mc: {}".format(result.value_at_initial, mc_res))

    print("iternum={}".format(solver.iterations))
    print("solver time={}".format(solver.solver_timer))


def get_prob01States(model, formulas):
    parameters = model.collect_probability_parameters()
    instantiator = stormpy.pars.PDtmcInstantiator(model)
    print(parameters)

    point = dict()
    for p in parameters:
        point[p] = stormpy.RationalRF(0.4)

        print(p)
    instantiated_model = instantiator.instantiate(point)
    assert instantiated_model.nr_states == model.nr_states
    assert not instantiated_model.has_parameters
    pathform = formulas[0].raw_formula.subformula
    assert type(pathform) == stormpy.logic.EventuallyFormula
    labelform = pathform.subformula
    labelprop = stormpy.core.Property("label-prop", labelform)
    phiStates = stormpy.BitVector(instantiated_model.nr_states, True)
    psiStates = stormpy.model_checking(instantiated_model, labelprop).get_truth_values()
    (prob0, prob1) = stormpy.compute_prob01_states(model, phiStates, psiStates)
    return prob0, prob1
    # (prob0A, prob1E) = stormpy.compute_prob01max_states(model, phiStates, psiSta


def find_rew0_states(property, model):
    formula = property.raw_formula
    assert type(formula) == stormpy.logic.RewardOperator
    path_formula = formula.subformula
    if type(path_formula) == stormpy.logic.EventuallyFormula:
        psi_formula = path_formula.subformula
    else:
        raise ValueError("Property type not supported")
    psi_result = stormpy.model_checking(model, psi_formula)
    psi_states = psi_result.get_truth_values()
    return psi_states
    # (prob0A, prob1E) = stormpy.compute_prob01max_states(model, phiStates, psiSta
if __name__ == '__main__':
    path = ""
    interval_path = ""
    formula_str = ""
    threshold = 320
    main(path, interval_path, formula_str, threshold)


