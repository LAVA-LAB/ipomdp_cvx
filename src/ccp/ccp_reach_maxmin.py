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
import ccp.interval_parser as interval_parser
#from util import interval_parser

#import value_iterator
import stormpy.pomdp
import time
import random




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



def build_robust_constraint(cons, coeff, denom1, pinit, succ, pVars, paraminit, param_id, paramVars):
    if coeff < 0:
        cons = cons + coeff * denom1 * (-0.5 * (pinit[succ]) ** 2 - pinit[succ] * (
                pVars[succ] - pinit[succ]) - 0.5 * (paraminit[param_id]) ** 2 -
                                        paraminit[
                                            param_id] * (
                                                paramVars[param_id] - paraminit[
                                            param_id]))
        cons = cons + coeff * denom1 * 0.5 * (pVars[succ] + paramVars[param_id]) * (
                pVars[succ] + paramVars[param_id])
        return cons
    else:
        cons = cons - coeff * denom1 * (-0.5 * (pinit[succ]) ** 2 - pinit[succ] * (
                pVars[succ] - pinit[succ]) - 0.5 * (paraminit[param_id]) ** 2 -
                                        paraminit[
                                            param_id] * (
                                                paramVars[param_id] - paraminit[
                                            param_id]))
        cons = cons - coeff * denom1 * 0.5 * (pVars[succ] - paramVars[param_id]) * (
                pVars[succ] - paramVars[param_id])
        return cons



class QcqpSolver():
    def __init__(self):
        self.solver_timer = 0.0
        self.encoding_timer = 0.0
        self.iterations = 0
        self.model_check_timer = 0.0

    def run(self, model, fsc_parameters, fsc_parameters_ids, interval_parameters, interval_parameters_ids, properties, prob0E, prob1A, threshold, direction, options, intervals, polyhedron_state_map, model_check=True):
        """
        Runs the QCQP procedure by a series of calls to gurobi.

        :param model: The model
        :type model: a stormpy dtmc/mdp
        :param parameters: The parameters occuring in the model
        :type parameters: a list of pycarl variables
        :param properties: The properties as an iterable over stormpy.properties
        :param threshold: The threshold
        :type threshold: float
        :param direction: Are we looking for a value below or above
        :type direction: a string, either "above" or "below"
        :param options: Further options with which the algorithm should run
        :return:
        """
        assert direction in ["above", "below"]
        if direction == "above":
            raise RuntimeError("Direction == above is currently not supported.")
        if not options.silent:
            print("Number of states: {}".format(model.nr_states))
            print("Number of transitions: {}".format(model.nr_transitions))
            print("Labels: {}".format(model.labeling.get_labels()))
            print(model.model_type)
            print("Number of states: {}".format(model.nr_states))


        numstate = model.nr_states
        # print(numstate)
        # Initializing some arrays for state, parameter and tau variables, and their values at previous iterations

        #paraminit = dict([[x.id, 0.5] for x in parameters if not x.name[0] == "I"])
        solution = dict()

        paraminit = {}
        for x in fsc_parameters:
            paraminit[x.id] = 0.5 # random.random()
            solution[x] = stormpy.RationalRF(paraminit[x.id])
        #for x in fsc_parameters:
        #    solution[x] = stormpy.RationalRF(paraminit[x.id])
        # print(solution)
        regiondict = dict()
        for x in interval_parameters:
            for interval in intervals:
                if interval.id == x.id:
                    # print("---\n"+str(type(x))+"\n---")
                    regiondict[x] = (
                        stormpy.RationalRF(interval.lowerbound), stormpy.RationalRF(interval.upperbound))
        region = stormpy.pars.ParameterRegion(regiondict)
        instantiator = stormpy.pars.PartialPDtmcInstantiator(model)
        instantiated_model = instantiator.instantiate(solution)

        env = stormpy.Environment()
        env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
        env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
        env.solver_environment.native_solver_environment.precision = stormpy.Rational('0.01')

        #env = stormpy.Environment()

        #env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
        #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration

        # old stuff i guess(?)
        # env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
        # env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration

        # env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
        # env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
        # env.solver_environment.native_solver_environment.precision = stormpy.Rational("0.01")
        start_check = time.time()

        region_checker = stormpy.pars.create_region_checker(env, instantiated_model,
                                                            properties[0].raw_formula,
                                                            allow_model_simplification=False)

        #print("region check")
        result = region_checker.get_bound_all_states(env, region, maximise=False)
        end_check = time.time()
        self.model_check_timer += (end_check - start_check)
        #print("model check time: " + str(model_check_timer))
        ansval = result.at(model.initial_states[0])

        results = []
        print("---------\n")
        print("{0},{1}".format(ansval, self.model_check_timer + self.encoding_timer + self.solver_timer))
        print("\n")
        results.append((ansval, self.model_check_timer + self.encoding_timer + self.solver_timer))


        #print("model checking init result: " + str(ansval))
        #for x in fsc_parameters:
        #    for y in interval_parameters:
        #        if x.id == y.id:
        #            print("EQUAL IDs")
        #print("EQUALITY CHECK DONE")

        #paraminit = dict([[x.id, 0.5]] for x in fsc_parameters)
        pinit = [threshold for _ in range(numstate)]
        for state in model.states:
            pinit[state.id] = (result.at(state))
        #cinit = [10.0 for _ in range(numstate)]

        # The penalty parameter for constraint violation
        mu = options.mu
        # Getting the initial state
        initstate = int(model.initial_states[0])
        for i in range(options.maxiter):
            self.iterations = i
            encoding_start = time.time()
            m = Model("qcp")
            m.setParam('OutputFlag', not options.silent)

            # Initializing some arrays for state, parameter and tau variables, and their values at previous iterations
            # Initializing gurobi variables for parameters,lb=lowerbound, ub=upperbound
            pVars = [m.addVar(lb=0, ub=1.0,name="p"+str(i)) for i in range(numstate)]
            #cVars = [m.addVar(lb=0, ub=1.0) for _ in range(numstate)]

            tau = [m.addVar(lb=0,name="t"+str(i)) for i in range(numstate)]
            #taucost = [m.addVar(lb=0) for _ in range(numstate)]

            #tt = m.addVar(lb=0.0, name="TT")

            #paramVars = dict([[x.id, m.addVar(lb=0)] for x in parameters if not x.name[0] == 'I'])

            paramVars = {}
            for x in fsc_parameters:
                paramVars[x.id] = m.addVar(lb=options.graph_epsilon,ub=1-options.graph_epsilon,name="param"+str(x.id))

            #print("ParamVars")
            #print(paramVars)
            #print("---")


            # A counter to check number of transitions
            numtrans = 0
            # List of constraints
            #constraints = []
            # Updates the model for gurobi
            m.update()
            #print(pVars)


            for state in model.states:
                #print(state.id)
                # print (prob1.get(state))
                # Cons=values constraints on the right hand side for a pdtmc

                # Find the polyhedrons at the current state
                current_polyhedrons = polyhedron_state_map[state.id]
                #for p in polyhedrons:
                #    #print(type(p.state))
                #    #print(type(state.id))
                #    if int(p.state) == int(state.id):
                #        current_polyhedrons.append(p)

                if len(current_polyhedrons) == 0:
                    # non robust constraint:
                    cons = 0
                    flag = 0
                    for action in state.actions:
                        for transition in action.transitions:

                            transition_value = transition.value()

                            numtrans = numtrans + 1

                            succ = int(transition.column)
                            # Value of transition
                            transition_value = transition.value()
                            # TODO check that the denominator is indeed constant?
                            # Denominator of transition
                            den = transition_value.denominator.constant_part()
                            denom1 = 1 / float(den)

                            # If the transition value is not constant
                            if not transition_value.is_constant():
                                num = transition_value.numerator.polynomial()

                                # Iterates over terms in numerators
                                for t in num:
                                    # If the transition term is a constant
                                    if t.is_constant():
                                        # Add just value of transition is the successor state is prob1 to the constraints
                                        if prob1A.get(succ):
                                            cons = cons + float(t.coeff) * denom1
                                        # Add nothing successor state is prob0
                                        elif prob0E.get(succ):
                                            pass
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            cons = cons + float(t.coeff) * denom1 * pVars[succ]

                                    # If the transition term is not a constant
                                    else:
                                        if t.tdeg > 1:
                                            raise RuntimeError("We expect the term to be a single variable")
                                        if t.monomial[0][0].name[0] == 'I':
                                            raise RuntimeError("Interval parameter in non-robust constraint")
                                        if t.monomial[0][0].id in interval_parameters_ids:
                                            print(t.monomial[0][0].id)
                                            print(t.monomial[0][0].name)
                                            print(interval_parameters_ids)
                                            print(state.id)
                                            print(current_polyhedrons)
                                            raise RuntimeError("Interval parameter in non-robust constraint")
                                        if t.monomial[0][0].id not in fsc_parameters_ids:
                                            print("param id = " + str(t.monomial[0][0].id))
                                            print("name = " + str(t.monomial[0][0].name))
                                            print(type(fsc_parameters))
                                            elem = fsc_parameters.pop()
                                            print(elem)
                                            print(type(elem))
                                            raise RuntimeError("Non fsc-parameter found in non-robust constraint")

                                        param_id = t.monomial[0][0].id
                                        coeff = float(t.coeff)

                                        # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                        if prob1A.get(succ):
                                            cons = cons + coeff * paramVars[param_id] * denom1
                                        # Add nothing successor state is prob0
                                        elif prob0E.get(succ):
                                            pass
                                        # Add the quadratic term to the constraints
                                        else:
                                            flag = 1

                                            # The bilinear terms are split into convex+concave terms, then the concave term is underapproximated by a affine term
                                            # First term in the addition is the affine term, second term is the convex term

                                            ######This portion changes for above compared to below
                                            if coeff < 0:
                                                cons = cons + coeff * denom1 * (-0.5 * (pinit[succ]) ** 2 - pinit[succ] * (pVars[succ] - pinit[succ]) -
                                                0.5 * (paraminit[param_id]) ** 2 - paraminit[param_id] * (paramVars[param_id] - paraminit[param_id]))
                                                cons = cons + coeff * denom1 * 0.5 * (pVars[succ] + paramVars[param_id]) * (pVars[succ] + paramVars[param_id])
                                            else:
                                                cons = cons - coeff * denom1 * (-0.5 * (pinit[succ]) ** 2 - pinit[succ] * (pVars[succ] - pinit[succ]) -
                                            0.5 * (paraminit[param_id]) ** 2 -paraminit[param_id] * (paramVars[param_id] - paraminit[param_id]))
                                                cons = cons - coeff * denom1 * 0.5 * (pVars[succ] - paramVars[param_id]) * (pVars[succ] - paramVars[param_id])

                            # If the value of transition is constant
                            else:
                                # Get the value of transition
                                constant_value = transition_value.constant_part()
                                # If successor state is prob1, just add the value of transition
                                if prob1A.get(succ):
                                    cons = cons + float(constant_value)
                                # If successor state is prob0, do nothing
                                elif prob0E.get(succ):
                                    pass
                                # Else, add transitionvalue*p_succ
                                else:
                                    cons = cons + float(constant_value) * pVars[succ]

                        # If the constraint is quadratic, add a penalty term to the constraints, otherwise dont add the term
                        ######These constraints are reversed for above
                        #print(cons,state.id)
                        if flag == 1:
                            m.addQConstr(pVars[state.id] <= cons + tau[state.id])
                        else:
                            m.addConstr(pVars[state.id] <= cons+ tau[state.id])
                        #print(state.id,pVars[state.id])



                else:
                    # robust constraint
                    # get all possible instantiations
                    combinations = interval_parser.calculate_vertex_combinations(current_polyhedrons)
                    #print(combinations,"combinations")
                    #print(len(combinations))
                    #print(combinations[0])
                    instantiation=combinations[0]
                    #print(instantiation,state.id)

                    # for each instantiation build a constraint
                    for instantiation in combinations:
                        for vertex in instantiation:
                            #print(instantiation,state.id)
                            #print(len(instantiation))
                            # for every instantiation, a constraint is added
                            cons = 0

                            # A flag for linear vs quadratic constraints
                            flag = 0

                            for action in state.actions:
                                for transition in action.transitions:

                                    transition_value = transition.value()
                                    numtrans = numtrans + 1
                                    succ = int(transition.column)
                                    # Value of transition
                                    transition_value = transition.value()
                                    # TODO check that the denominator is indeed constant?
                                    # Denominator of transition
                                    den = transition_value.denominator.constant_part()
                                    denom1 = 1 / float(den)

                                    # If the transition value is not constant
                                    if not transition_value.is_constant():
                                        num = transition_value.numerator.polynomial()

                                        # Iterates over terms in numerators
                                        for t in num:
                                            # If the transition term is a constant
                                            if t.is_constant():
                                                # Add just value of transition is the successor state is prob1 to the constraints
                                                if prob1A.get(succ):
                                                    cons = cons + float(t.coeff) * denom1
                                                # Add nothing successor state is prob0
                                                elif prob0E.get(succ):
                                                    pass
                                                # Else add transitionvalue*p_succ to the constraint
                                                else:
                                                    cons = cons + float(t.coeff) * denom1 * pVars[succ]

                                            # If the transition term is not a constant
                                            else:
                                                #if t.tdeg > 1:
                                                #    raise RuntimeError("We expect the term to be a single variable")
                                                #print(t)
                                                if t.tdeg == 1:
                                                    # term of degree 1, check if it is an interval or a real parameter
                                                    param_name = t.monomial[0][0].name
                                                    param_id = t.monomial[0][0].id

                                                    #if param_name[0] == 'I':
                                                    if param_id in interval_parameters_ids:

                                                        #interval_name = param_name
                                                        interval_id = param_id
                                                        interval_value = -1
                                                        #print(instantiation)
                                                            #print(len(instantiation))
                                                        for vertex_value in vertex.vertex_values:
                                                                #print(type(vertex_value.id))
                                                                #print(type(interval_id))
                                                            if vertex_value.id == interval_id:
                                                                interval_value = vertex_value.value
                                                                #print(interval_value,interval_id,t.monomial[0][0].name)

                                                                coeff = float(t.coeff) * float(interval_value)

                                                                break
                                                    if prob1A.get(succ):

                                                    # no real parameter, only interval coeff.
                                                        cons = cons + coeff * denom1
                                                    elif prob0E.get(succ):
                                                        pass
                                                    else:
                                                        cons = cons + coeff * denom1* pVars[succ]
                                                    if interval_value == -1:
                                                        print(interval_value,state.id)
                                                        raise RuntimeError("Interval not part of the instantiation")
                                                    else:
                                                        param_id = t.monomial[0][0].id
                                                        coeff = float(t.coeff)
                                                elif t.tdeg >=2:
                                                    raise RuntimeError("Does not make sense for simple POMDPs, the degree cannot be 2")
                                                    par1 = t.monomial[0][0].id
                                                    par2 = t.monomial[1][0].id
                                                    if par1 in interval_parameters_ids:
                                                        #print(par1)
                                                        param_id = t.monomial[1][0].id
                                                        interval_id = par1
                                                    elif par2 in interval_parameters_ids:
                                                        #print(par2)
                                                        param_id = t.monomial[0][0].id
                                                        interval_id = par2
                                                    else:
                                                        raise RuntimeError("Multiple non-interval parameters in the same monomial")

                                                    # find value of the interval given instantiation
                                                    interval_value = -1
                                                    for vertex_value in vertex.vertex_values:
                                                        if vertex_value.id == interval_id:
                                                            interval_value = vertex_value.value
                                                            coeff = float(t.coeff)*float(interval_value)
                                                            break
                                                    if interval_value == -1:
                                                        #print("interval name: " + str(interval_name))
                                                        raise RuntimeError("Interval not part of the instantiation")

                                                else:
                                                    raise RuntimeError("We expect the term to be a single variable")

                                                # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                                if param_id == None:
                                                    if prob1A.get(succ):

                                                    # no real parameter, only interval coeff.
                                                        cons = cons + coeff * denom1
                                                    elif prob0E.get(succ):
                                                        pass
                                                    else:
                                                        cons = cons + coeff * denom1* pVars[succ]


                                    # If the value of transition is constant
                                    else:
                                        # Get the value of transition
                                        constant_value = transition_value.constant_part()
                                        # If successor state is prob1, just add the value of transition
                                        if prob1A.get(succ):
                                            cons = cons + float(constant_value)
                                        # If successor state is prob0, do nothing
                                        elif prob0E.get(succ):
                                            pass
                                        # Else, add transitionvalue*p_succ
                                        else:
                                            cons = cons + float(constant_value) * pVars[succ]

                                # If the constraint is quadratic, add a penalty term to the constraints, otherwise dont add the term
                                ######These constraints are reversed for above
                                #print(cons,state.id)
                                if flag == 1:
                                    m.addQConstr(pVars[state.id] <= cons + tau[state.id])
                                else:
                                    m.addConstr(pVars[state.id] <= cons+ tau[state.id])
                            #print(state.id,pVars[state.id])

            ##Adds constraints for prob1 and prob0 state
        #    for state in range(numstate):
        #        if prob1A.get(state):
        #            m.addConstr(pVars[state] == 1)
        #        elif prob0E.get(state):
        #            m.addConstr(pVars[state] == 0)

            # Constraint for initial state
            ######This constraint is reversed for above

            m.addConstr(pVars[initstate] >= threshold)
            objective = 0.0
            # Adding terms to the objective
            for state in range(numstate):
                # This constraints minimizes the max of violation
                #m.addConstr(tau[state] <= tt)
                # This minimizes sum of violation, mu is the parameter that punishes violation
                objective = objective + tau[state]
                objective = objective - pVars[state]/mu

            # Maximize the probability of initial state
      #      objective = objective - pVars[initstate]
            self.encoding_timer += (time.time() - encoding_start)

            start3 = time.time()
            # Solves the problem
            m.setObjective(objective, GRB.MINIMIZE)


            # write model to a file for debugging
            #m.write("last_model.mps")

            #print('Solving...')
            m.optimize()

            t3 = time.time()
            self.solver_timer += (t3 - start3)

            #print("Solver time :" + str(t3 - start3))
            # Prints the maximum violation
            maxx = 0
            for state in range(numstate):
                val = tau[state].x
                if val > maxx:
                    maxx = val

            #if not options.silent:
            #    print("Max vio :", maxx)
            #    print("p =", pVars[initstate].x)

            # Breaks if the violation is small and prints number of iterations and total time

            # Updates the probability values for next iteration


            # Updares the parameter values for next iteration
            #print(paramVars)
            for param_id, param_var in paramVars.items():
                if not isinstance(param_var, int):
                    if abs(param_var.x) > options.graph_epsilon:
                        #  print pVar
                        paraminit[param_id] = param_var.x
                    else:
                        paraminit[param_id] = options.graph_epsilon

            #parameters = model.collect_probability_parameters()
            #parameter_values = dict([[id, param_var.x] for id, param_var in paramVars.items()])
            #parameter_assignments = dict([[x.name, parameter_values[x.id]] for x in fsc_parameters])
            parameter_names = dict([[x.id, x.name] for x in fsc_parameters])

            param_values = {}
            for p in fsc_parameters_ids:
                param_values[p] = paramVars[p].x
            #print(param_values)
            #param_values = dict([[id, param_var.x] for id, param_var in paramVars.items()])


            # model checking for early termination:

            maxx = 0
            trust_region = 2.5
            bestval = pinit[initstate]
            model_check_timer = 0.0

            if model_check:
                # print(solution)
                option = "max"  # min or max for value iterator
                solution = dict()
                for x in fsc_parameters:
                    solution[x] = stormpy.RationalRF(param_values[x.id])
                # print(solution)
                regiondict = dict()
                for x in interval_parameters:
                    for interval in intervals:
                        if interval.id == x.id:
                            # print("---\n"+str(type(x))+"\n---")
                            regiondict[x] = (
                                stormpy.RationalRF(interval.lowerbound), stormpy.RationalRF(interval.upperbound))
                region = stormpy.pars.ParameterRegion(regiondict)
                instantiator = stormpy.pars.PartialPDtmcInstantiator(model)
                instantiated_model = instantiator.instantiate(solution)
                # assert interval_parameters == instantiated_model.collect_probability_parameters()

                env = stormpy.Environment()
                env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
                env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
                env.solver_environment.native_solver_environment.precision = stormpy.Rational('0.01')


                # settings 1
                #env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
                #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.optimistic_value_iteration
                #env.solver_environment.minmax_solver_environment.precision = stormpy.Rational("0.01")

                # settings 2
                #env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
                #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration
                #env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
                #env.solver_environment.native_solver_environment.precision = stormpy.Rational("0.01")

                # settings 3
                # env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
                # env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration

                # old stuff i guess(?)
                # env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
                # env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration

                # env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
                # env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
                # env.solver_environment.native_solver_environment.precision = stormpy.Rational("0.01")
                start_check = time.time()

                region_checker = stormpy.pars.create_region_checker(env, instantiated_model,
                                                                    properties[0].raw_formula,
                                                                    allow_model_simplification=False)

                #print("region check")
                result = region_checker.get_bound_all_states(env, region, maximise=False)
                end_check = time.time()
                self.model_check_timer += (end_check - start_check)
                #print("model check time: " + str(model_check_timer))

                # paraminit = dict([[x.id, float(stormpy.RationalRF(solution[x]))] for x in self._parameters])
                for x in fsc_parameters:
                    paraminit[x.id] = float(stormpy.RationalRF(solution[x]))

                ansval = result.at(model.initial_states[0])
                totaltime = self.solver_timer + self.encoding_timer + self.model_check_timer
                print("{0}, {1}".format(ansval, totaltime))
                results.append((ansval, totaltime))

                if abs(maxx) < options.graph_epsilon:
                    print("termination due to slack variables being epsilon close to zero at iteration {0}: ".format(str(i)))
                    return results, "slack variables epsilon close to zero"

                if ansval > threshold:
                    print(
                        "Early termination due to positive model checking result at iteration {0}: ".format(str(i)))
                    return results, "policy evaluation passed threshold"
                elif totaltime > options.timeout:
                    print(
                        "termination due to timeout: {0}".format(str(totaltime)))
                    return results, "timeout"
                else:
                    # print(pinit)
                    for state in model.states:
                        pinit[state.id] = (result.at(state))
                    # Updares the parameter values for next iteration
                    for param_id, param_var in paramVars.items():
                        if not isinstance(param_var, int):
                            if abs(param_var.x) > options.graph_epsilon:
                                #  print pVar
                                paraminit[param_id] = param_var.x
                            else:
                                paraminit[param_id] = options.graph_epsilon
                    self.solver_params = solution

                m.update()

                #print("bestval:", bestval)
            maxp=0
            for state in range(numstate):
                 #if abs(cVars[state].x) > options.graph_epsilon:
                #pinit[state] = mc_res1.at(state)
                if 1>pinit[state]>maxp:
                    maxp=pinit[state]
            # Updates penalty parameter
            mu = min(mu * (maxp+1), 1e10)


        return results, "max iterations"





def main_prism(path, interval_path, formula_str, threshold, memval=1,timeout=1800, maxiter=200, evaluation_set = []):
    prism_program = stormpy.parse_prism_program(path)

    # formula_str = "P=? [!\"bad\" U \"goal\"]"
    # formula_str = "P=? [F \"goal\"]"
    # interval_path="collision_partial_obs_2d_upd_hobs_20_small.intervals"
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    properties = stormpy.parse_properties_for_prism_program(formula_str, prism_program)
    # construct the pPOMDP
    pomdp = stormpy.build_parametric_model(prism_program, properties)
    pomdp_parameters = pomdp.collect_probability_parameters()
    return main(pomdp, interval_path, formula_str, threshold, memval, path, timeout, maxiter, evaluation_set)


def main_drn(path, interval_path, formula_str, threshold, memval=1,timeout=1800, maxiter=200, evaluation_set = []):
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    pomdp = stormpy.build_parametric_model_from_drn(path, opts)
    return main(pomdp, interval_path, formula_str, threshold, memval, path, timeout, maxiter, evaluation_set)



def main(pomdp, interval_path, formula_str, threshold, memval=1, path="", timeout=1800, maxiter=200, evaluation_set = []):
    model_info = dict()
    model_info["model name"] = path
    model_info["objective"] = formula_str

    t0 = time.time()

    pomdp_parameters = pomdp.collect_probability_parameters()
    # stormpy.export_to_drn(pomdp, "pomdp_ex")
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



    for p in polyhedrons:
        p.compute_vertices()


    properties = stormpy.parse_properties(formula_str)

    prob0E, prob1A = stormpy.prob01max_states(pmc, properties[0].raw_formula.subformula)

    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=0.05, maxiter=maxiter, graph_epsilon=1e-2, silent=True, timeout=timeout)


    solver = QcqpSolver()


    fsc_parameters_ids = set()
    for x in fsc_parameters:
        fsc_parameters_ids.add(x.id)
    interval_parameters_ids = set()
    for x in pomdp_parameters:
        interval_parameters_ids.add(x.id)

    polyhedron_state_map = dict()
    for id in range(pmc.nr_states):
        current_polyhedrons = []
        for p in polyhedrons:
            if int(p.state) == int(id):
                current_polyhedrons.append(p)
        polyhedron_state_map[id] = current_polyhedrons

    solver_results, solver_exit = solver.run(pmc, fsc_parameters, fsc_parameters_ids, pomdp_parameters, interval_parameters_ids, properties, prob0E, prob1A, threshold, direction, options, intervals, polyhedrons, polyhedron_state_map, True)





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

    tend = time.time()

    print("Final result: ", eval_results[interval_path])
    print("Total time: ", str(tend - t0))

    model_info["total computation time"] = str(tend - t0)
    model_info["final result"] = eval_results[interval_path]

    return model_info, solver_results, eval_results






    """
    #compute the policy against the robust interval
    #interval_path="collision_partial_obs_2d_upd_hobs_20_small.intervals"
    #intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)
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
    ansval=result.at(pmc.initial_states[0])
    tend = time.time()

    print("\n\n----------")
    print("Final result: ", ansval)
    print("Total time: ", str(tend - t0))

    model_info["total computation time"] = str(tend - t0)
    model_info["final result"] = ansval

    return model_info, solver_results
    """




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

















