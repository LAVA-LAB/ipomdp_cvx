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
#import value_iterator
import stormpy.pomdp
import time
import random



def input_files():
    model_path = "collision_partial_obs_2d_upd_hobs_20_small_2.drn"
    interval_path = "collision_partial_obs_2d_upd_hobs_20_big.intervals"
    formula_str = "P=? [F \"goal\"]"
    threshold = 0.75

    return model_path, interval_path, formula_str, threshold


class QcqpOptions():
    def __init__(self, mu, maxiter, graph_epsilon, silent):
        self.mu = mu
        self.maxiter = maxiter
        self.graph_epsilon = graph_epsilon
        self.silent = silent


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

    def run(self, model, fsc_parameters, fsc_parameters_ids, interval_parameters, interval_parameters_ids, properties, prob0E, prob1A, threshold, direction, options, intervals, polyhedrons, polyhedron_state_map, model_check):
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
        for x in fsc_parameters:
            solution[x] = stormpy.RationalRF(paraminit[x.id])
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
        env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration

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

        print("region check")
        result = region_checker.get_bound_all_states(env, region, maximise=False)
        end_check = time.time()
        model_check_timer = end_check - start_check
        print("model check time: " + str(model_check_timer))
        ansval = result.at(model.initial_states[0])
        print("model checking init result: " + str(ansval))
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

                                                #else:
                                                    #raise RuntimeError("On a robust constraint, no paramVars should happen")

                                                    #if prob1A.get(succ):
                                                    #    cons = cons + coeff * paramVars[param_id] * denom1


                                                    # Add nothing successor state is prob0
                                                    # elif prob0E.get(succ):
                                                    #     pass
                                                    # # Add the quadratic term to the constraints
                                                    # else:
                                                    #     flag = 1
                                                    #
                                                    #     # The bilinear terms are split into convex+concave terms, then the concave term is underapproximated by a affine term
                                                    #     # First term in the addition is the affine term, second term is the convex term
                                                    #
                                                    #     ######This portion changes for above compared to below
                                                    #     if coeff < 0:
                                                    #         cons = cons + coeff * denom1 * (-0.5 * (pinit[succ]) ** 2 - pinit[succ] * (
                                                    #             pVars[succ] - pinit[succ]) - 0.5 * (paraminit[param_id]) ** 2 -
                                                    #                                         paraminit[
                                                    #                                             param_id] * (
                                                    #                                             paramVars[param_id] - paraminit[
                                                    #                                                 param_id]))
                                                    #         cons = cons + coeff * denom1 * 0.5 * (pVars[succ] + paramVars[param_id]) * (
                                                    #             pVars[succ] + paramVars[param_id])
                                                    #     else:
                                                    #         cons = cons - coeff * denom1 * (-0.5 * (pinit[succ]) ** 2 - pinit[succ] * (
                                                    #             pVars[succ] - pinit[succ]) - 0.5 * (paraminit[param_id]) ** 2 -
                                                    #                                         paraminit[
                                                    #                                             param_id] * (
                                                    #                                             paramVars[param_id] - paraminit[
                                                    #                                                 param_id]))
                                                    #         cons = cons - coeff * denom1 * 0.5 * (pVars[succ] - paramVars[param_id]) * (
                                                    #             pVars[succ] - paramVars[param_id])

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

            """ REMOVE TRIVIAL CONSTRAINTS FOR SIMPLE POMDPs
            for state in model.states:
                # Find the polyhedrons at the current state
                current_polyhedrons = []
                for p in polyhedrons:
                    if int(p.state) == int(state.id):
                        current_polyhedrons.append(p)

                if len(current_polyhedrons) == 0:
                    # non robust constraint
                    # cons = 0
                    for action in state.actions:

                        for transition in action.transitions:
                            cons1 = 0.0
                            transition_value = transition.value()
                            den = transition_value.denominator.constant_part()
                            # If the transition value is not a constant
                            if not transition_value.is_constant():
                                num = transition_value.numerator.polynomial()

                                for t in num:
                                    # Adds the value of the transition
                                    if t.is_constant():
                                        cons1 = cons1 + float(t.coeff) / float(den)
                                    # Adds the value of the transition
                                    else:
                                        if t.tdeg > 1:
                                            raise RuntimeError("We expect the term to be a single variable")
                                        coeff = float(t.coeff)
                                        #print(t.monomial[0][0].name)
                                        param_id = t.monomial[0][0].id
                                        cons1 = cons1 + coeff * paramVars[param_id] / float(den)
                            # If the transition has parameters, constrain each transition to between 0 and 1
                            if not isinstance(cons1, float):
                                # print(cons1)
                                m.addConstr(cons1 >= 0 + options.graph_epsilon)
                                m.addConstr(cons1 <= 1 - options.graph_epsilon)


                else:
                    # robust constraint

                    # get all possible instantiations
                    combinations = interval_parser.calculate_vertex_combinations(current_polyhedrons)
                    # for each instantiation build a constraint
                    for instantiation in combinations:

                        #cons = 0
                        for action in state.actions:

                            for transition in action.transitions:
                                cons1 = 0.0
                                transition_value = transition.value()
                                den = transition_value.denominator.constant_part()
                                # If the transition value is not a constant
                                if not transition_value.is_constant():
                                    num = transition_value.numerator.polynomial()

                                    for t in num:
                                        # Adds the value of the transition
                                        if t.is_constant():
                                            cons1 = cons1 + float(t.coeff) / float(den)
                                        # Adds the value of the transition
                                        else:

                                            # if t.tdeg > 1:
                                            #    raise RuntimeError("We expect the term to be a single variable")
                                            if t.tdeg == 1:
                                                param_name = t.monomial[0][0].name
                                                param_id = t.monomial[0][0].id
                                                if param_id in interval_parameters_ids:

                                                    interval_id = param_id
                                                    interval_value = -1
                                                    for vertex in instantiation:
                                                        for vertex_value in vertex.vertex_values:
                                                            if vertex_value.id == interval_id:
                                                                interval_value = vertex_value.value
                                                                coeff = float(t.coeff) * float(interval_value)
                                                                break
                                                    if interval_value == -1:
                                                        raise RuntimeError("Interval not part of the instantiation")
                                                    param_id = None
                                                else:
                                                    param_id = t.monomial[0][0].id
                                                    coeff = float(t.coeff)
                                            elif t.tdeg == 2:
                                                par1 = t.monomial[0][0].id
                                                par2 = t.monomial[1][0].id
                                                if par1 in interval_parameters_ids:
                                                    param_id = t.monomial[1][0].id
                                                    interval_id = par1
                                                elif par2 in interval_parameters:
                                                    param_id = t.monomial[0][0].id
                                                    interval_id = par2
                                                else:
                                                    raise RuntimeError("Multiple non-interval parameters in the same monomial")

                                                # find value of the interval given instantiation
                                                interval_value = -1
                                                for vertex in instantiation:
                                                    for vertex_value in vertex.vertex_values:
                                                        if vertex_value.id == interval_id:
                                                            interval_value = vertex_value.value
                                                            coeff = float(t.coeff) * interval_value
                                                            break
                                                if interval_value == -1:
                                                    print(interval_name)
                                                    raise RuntimeError("Interval not part of the instantiation")

                                            else:
                                                raise RuntimeError("We expect the term to be a single variable")

                                            if param_id == None:
                                                cons1 = cons1 + coeff / float(den)
                                            else:
                                                cons1 = cons1 + coeff * paramVars[param_id] / float(den)
                                # If the transition has parameters, constrain each transition to between 0 and 1
                                if not isinstance(cons1, float):
                                    # print(cons1)
                                    m.addConstr(cons1 >= 0 + options.graph_epsilon)
                                    m.addConstr(cons1 <= 1 - options.graph_epsilon)
            """

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

            print('Solving...')
            m.optimize()

            t3 = time.time()
            self.solver_timer += (t3 - start3)

            print("Solver time :" + str(t3 - start3))
            # Prints the maximum violation
            maxx = 0
            for state in range(numstate):
                val = tau[state].x
                if val > maxx:
                    maxx = val

            if not options.silent:
                print("Max vio :", maxx)
                print("p =", pVars[initstate].x)

            # Breaks if the violation is small and prints number of iterations and total time
            # if abs(maxx) < options.graph_epsilon:
            #     param_values = dict([[id, param_var.x] for id, param_var in paramVars.items()])
            #     return QcqpResult(pVars[initstate].x, param_values)
            # Updates the probability values for next iteration
            for state in range(numstate):
                pass
                #if abs(pVars[state].x) > options.graph_epsilon:
                #    pinit[state] = pVars[state].x
                #else:
                #    pinit[state] = pVars[state].x

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

                # settings 1
                #env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
                #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.optimistic_value_iteration
                #env.solver_environment.minmax_solver_environment.precision = stormpy.Rational("0.01")

                # settings 2
                env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.native)
                env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration
                env.solver_environment.native_solver_environment.method = stormpy.NativeLinearEquationSolverMethod.optimistic_value_iteration
                env.solver_environment.native_solver_environment.precision = stormpy.Rational("0.01")

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

                print("region check")
                result = region_checker.get_bound_all_states(env, region, maximise=False)
                end_check = time.time()
                model_check_timer = end_check - start_check
                print("model check time: " + str(model_check_timer))

                # paraminit = dict([[x.id, float(stormpy.RationalRF(solution[x]))] for x in self._parameters])
                for x in fsc_parameters:
                    paraminit[x.id] = float(stormpy.RationalRF(solution[x]))

                ansval = result.at(model.initial_states[0])
                print("model checking result: " + str(ansval))
                # terminate if the result is better than threshold
                if ansval > threshold:
                    # TODO adjust result
                    print(
                        "Early termination due to positive model checking result at iteration {0}: ".format(str(i)))
                    print("p[init] = " + str(ansval))
                    print("SCP parameter values: ")
                    # for id, param_var in paramVars.items():
                    #    print(str(parameter_names[id]) + "  :  " + str(param_var.x))
                    # print("total model check time:", self.model_check_timer)
                    solver_params = solution

                    return QcqpResult(pVars[initstate].x, solver_params)
                # enlarge trust region and update probability and parameter values
                # if the best found solution is better

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

                print("bestval:", bestval)
            maxp=0
            for state in range(numstate):
                 #if abs(cVars[state].x) > options.graph_epsilon:
                #pinit[state] = mc_res1.at(state)
                if 1>pinit[state]>maxp:
                    maxp=pinit[state]
            #mu = mu * (maxp+1)
            # Updates penalty parameter
            mu = mu * (maxp+1)
            if mu>1e10:
                mu=1e10
            print("mu:",mu)

            """ OLD MODEL CHECKER
            if model_check:
                print("Model checking")
                solution = param_values
                option = "max" # min or max for value iterator
                ans = value_iterator.value_iteration_iterator(model, fsc_parameters, solution, interval_parameters_ids, prob0E, prob1A, polyhedrons, option)
                print("model checking ans:")
                print(ans)
                if ans < threshold:
                    # TODO adjust result
                    print("Early termination due to positive model checking result at iteration {0}: ".format(str(i)))
                    print("p[init] = " + str(ans) + " < " + str(threshold))
                    print("QCQP parameter values: ")
                    for id, param_var in paramVars.items():
                        print(str(parameter_names[id]) + "  :  " + str(param_var.x))
                    return QcqpResult(pVars[initstate].x, param_values)



            # model checking result failed:

            # bounds reached, end
            if pVars[initstate].x < threshold or i >= options.maxiter-1:
                print("p[init] = " + str(pVars[initstate].x))
                print("parameter values: ")
                #for id, param_var in paramVars.items():
                #    print(str(parameter_names[id]) + "  :  " + str(param_var.x))
                return QcqpResult(pVars[initstate].x, param_values)
            # remove storm model checking?

            #instantiator = stormpy.pars.PDtmcInstantiator(model)

            # # Check distance to result by storm (notice that also the point that storm checks slightly differs)
            #rational_parameter_assignments = dict([[x, stormpy.RationalRF(val)] for x, val in parameter_assignments.items()])
          #  print(rational_parameter_assignments)
            #print(type(rational_parameter_assignments))
            #instantiated_model = instantiator.instantiate(rational_parameter_assignments)
            #mc_res = stormpy.model_checking(instantiated_model, properties[0]).at(model.initial_states[0])
            #mc_res1 = stormpy.model_checking(instantiated_model, properties[0])
            #print("Mc: {}".format( mc_res))
            #if  mc_res>=threshold:
            #    param_values = dict([[id, param_var.x] for id, param_var in paramVars.items()])
            #    return QcqpResult(cVars[initstate].x, param_values)

            #maxp=0
            #for state in range(numstate):
                 #if abs(cVars[state].x) > options.graph_epsilon:
            #    pinit[state] = mc_res1.at(state)
            #    if 1>pinit[state]>maxp:
            #        maxp=pinit[state]
            #mu = mu * (maxp+1)
            ## Updates penalty parameter
            #mu = mu * (maxp+1)
            #if mu>1e10:
            #    mu=1e10


            """




# TODO make this also robust?
class QcqpRewSolver():
    def __init__(self):
        self.solver_timer = 0.0
        self.encoding_timer = 0.0
        self.iterations = 0

    def run(self, reward_model_name, model_rew, parameters_rew, rew0, threshold, direction, options, model_check):
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
            print("Number of states: {}".format(model_rew.nr_states))
            print("Number of transitions: {}".format(model_rew.nr_transitions))
            print("Labels: {}".format(model_rew.labeling.get_labels()))
            print(model_rew.model_type)
            print("Number of states: {}".format(model_rew.nr_states))


        numstate = model_rew.nr_states
        # print(numstate)
        # Initializing some arrays for state, parameter and tau variables, and their values at previous iterations
        paraminit = dict([[x.id, 0.5] for x in parameters_rew])

     #   pinit = [0.5 for _ in range(numstate)]
        cinit = [threshold for _ in range(numstate)]

        # The penalty parameter for constraint violation
        mu = options.mu
        # Getting the initial state
        initstate = int(model_rew.initial_states[0])
        for i in range(options.maxiter):
            self.iterations = i
            encoding_start = time.time()
            m = Model("qcp")
            m.setParam('OutputFlag', not options.silent)

            # Initializing some arrays for state, parameter and tau variables, and their values at previous iterations
            # Initializing gurobi variables for parameters,lb=lowerbound, ub=upperbound
        #    pVars = [m.addVar(lb=0, ub=1.0) for _ in range(numstate)]
            cVars = [m.addVar(lb=0) for _ in range(numstate)]

            tau = [m.addVar(lb=0) for _ in range(numstate)]
         #   taucost = [m.addVar(lb=0) for _ in range(numstate)]

            #tt = m.addVar(lb=0.0, name="TT")

            paramVars = dict([[x.id, m.addVar(lb=0)] for x in parameters_rew])
            #print(parameters_rew)

            # A counter to check number of transitions
            numtrans = 0
            # List of constraints
            #constraints = []
            # Updates the model for gurobi
            m.update()


            for state in model_rew.states:
                #assert model_rew.reward_models[reward_model_name].has_state_action_rewards
              #  print(reward_at_state)
                # print (prob1.get(state))
                # Cons=values constraints on the right hand side for a pdtmc
                cons = 0
                # A flag for linear vs quadratic constraints
                flag = 0

                for action in state.actions:
                    reward_at_state = model_rew.reward_models[reward_model_name].state_rewards[int(state.id)]
         #           print(reward_at_state)
                    if not reward_at_state.is_constant():
        #                print(reward_at_state.numerator.polynomial())
                        for term in reward_at_state.numerator.polynomial():
          #                  print(term.monomial[0][0].id)
                            if not term.is_constant():
                                param_id=term.monomial[0][0].id

                                cons=cons+paramVars[param_id]*float(term.coeff)
                            else:
                                cons=cons+float(term.coeff)
                    else:
                        cons=cons+float(reward_at_state.constant_part())

                    for transition in action.transitions:
                     #   print("From state {}, with probability {}, go to state {}".format(state, transition.value(), transition.column))

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
                                    if rew0.get(succ):
                                        pass
                                    # Add nothing successor state is prob0
                                    # elif prob0E.get(succ):
                                    #     pass
                                    # Else add transitionvalue*p_succ to the constraint
                                    else:
                                        cons = cons + float(t.coeff) * denom1 * cVars[succ]

                                # If the transition term is not a constant
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    param_id = t.monomial[0][0].id

#                                    print(param_id)
                                    coeff = float(t.coeff)
                                    # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                    if rew0.get(succ):
                                        pass
                                    # Add nothing successor state is prob0
                                    # elif prob0E.get(succ):
                                    #     pass
                                    # Add the quadratic term to the constraints
                                    else:
                                        flag = 1

                                        # The bilinear terms are split into convex+concave terms, then the concave term is underapproximated by a affine term
                                        # First term in the addition is the affine term, second term is the convex term

                                        ######This portion changes for above compared to below
                                        if coeff < 0:
                                            cons = cons + coeff * denom1 * (-0.5 * (cinit[succ]) ** 2 - cinit[succ] * (cVars[succ] - cinit[succ]) - 0.5 * (paraminit[param_id]) ** 2 -paraminit[param_id] * (paramVars[param_id] - paraminit[param_id]))
                                            cons = cons + coeff * denom1 * 0.5 * (cVars[succ] + paramVars[param_id]) * (cVars[succ] + paramVars[param_id])
                                        else:
                                            cons = cons - coeff * denom1 * (-0.5 * (cinit[succ]) ** 2 - cinit[succ] * (cVars[succ] - cinit[succ]) - 0.5 * (paraminit[param_id]) ** 2 -paraminit[param_id] * (paramVars[param_id] - paraminit[param_id]))
                                            cons = cons - coeff * denom1 * 0.5 * (cVars[succ] - paramVars[param_id]) * (cVars[succ] - paramVars[param_id])

                        # If the value of transition is constant
                        else:
                            # Get the value of transition
                            constant_value = transition_value.constant_part()
                            # If successor state is prob1, just add the value of transition
                            if rew0.get(succ):
                                pass
                            # If successor state is prob0, do nothing
                            # elif prob0E.get(succ):
                            #     pass
                            # Else, add transitionvalue*p_succ
                            else:
                                cons = cons + float(constant_value) * cVars[succ]

                    # If the constraint is quadratic, add a penalty term to the constraints, otherwise dont add the term

                    ######These constraints are reversed for above
                    if flag == 1:
                        m.addQConstr(cVars[state.id] <= cons + tau[state.id])
                    else:
                        m.addConstr(cVars[state.id] <= cons)
                    #print(state.id)
            for state in model_rew.states:

                #cons = 0
                for action in state.actions:

                    for transition in action.transitions:
                        cons1 = 0.0
                        transition_value = transition.value()
                        den = transition_value.denominator.constant_part()
                        # If the transition value is not a constant
                        if not transition_value.is_constant():
                            num = transition_value.numerator.polynomial()

                            for t in num:
                                # Adds the value of the transition
                                if t.is_constant():
                                    cons1 = cons1 + float(t.coeff) / float(den)
                                # Adds the value of the transition
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    param_id = t.monomial[0][0].id
                                    coeff = float(t.coeff)

                                    cons1 = cons1 + coeff * paramVars[param_id] / float(den)
                        # If the transition has parameters, constrain each transition to between 0 and 1
                        if not isinstance(cons1, float):
                            # print(cons1)
                            m.addConstr(cons1 >= 0 + options.graph_epsilon)
                            m.addConstr(cons1 <= 1 - options.graph_epsilon)
            ##Adds constraints for prob1 and prob0 state
        #    for state in range(numstate):
        #        if prob1A.get(state):
        #            m.addConstr(pVars[state] == 1)
        #        elif prob0E.get(state):
        #            m.addConstr(pVars[state] == 0)

            # Constraint for initial state
            ######This constraint is reversed for above

            m.addConstr(cVars[initstate] >= threshold)
            objective = 0.0
            # Adding terms to the objective
            for state in range(numstate):
                # This constraints minimizes the max of violation
                #m.addConstr(tau[state] <= tt)
                # This minimizes sum of violation, mu is the parameter that punishes violation
                objective = objective + tau[state]
                objective = objective - cVars[state]/mu


            self.encoding_timer += (time.time() - encoding_start)

            start3 = time.time()
            # Solves the problem
            m.setObjective(objective, GRB.MINIMIZE)
            print('Solving...')
            m.optimize()

            t3 = time.time()
            self.solver_timer += (t3 - start3)

            print("Solver time :" + str(t3 - start3))
            # Prints the maximum violation
            maxx = 0
            for state in range(numstate):
                val = tau[state].x
                if val > maxx:
                    maxx = val

            if not options.silent:
                print("Max vio :", maxx)
                print("p =", cVars[initstate].x)

            # Breaks if the violation is small and prints number of iterations and total time
            if abs(maxx) < options.graph_epsilon:
                param_values = dict([[id, param_var.x] for id, param_var in paramVars.items()])
                return QcqpResult(cVars[initstate].x, param_values)
            # Updates the probability values for next iteration
            for state in range(numstate):
                if abs(cVars[state].x) > options.graph_epsilon:
                    cinit[state] = cVars[state].x
                else:
                    cinit[state] = cVars[state].x
                #print(cinit[state],state)

            # Updares the parameter values for next iteration
            for param_id, param_var in paramVars.items():
                if not isinstance(param_var, int):
                    if abs(param_var.x) > options.graph_epsilon:
                        #  print pVar
                        paraminit[param_id] = param_var.x
                    else:
                        paraminit[param_id] = options.graph_epsilon
            # Updates penalty parameter
            mu = mu * 10.0
            if mu>1e8:
                mu=1e8



def main(model_path, interval_path, formula_str, threshold):
    t0 = time.time()

    #model_path = "collision_partial_obs_2d_upd_hobs_20_small_2.drn"
    #interval_path = "collision_partial_obs_2d_upd_hobs_20_small.intervals"
    #formula_str = "P=? [F \"goal\"]"
    #model_path, interval_path, formula_str, threshold = input_files()

    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    pomdp_drn_model = stormpy.build_parametric_model_from_drn(model_path, opts)

    interval_parameters = pomdp_drn_model.collect_probability_parameters()

    #print(interval_parameters)

    pomdp = stormpy.pomdp.make_canonic(pomdp_drn_model)


    # construct the memory for the FSC
    # in this case, a selective counter with two states
    memory_builder = stormpy.pomdp.PomdpMemoryBuilder()
    #number of memory states
    memval=1
    memory = memory_builder.build(stormpy.pomdp.PomdpMemoryPattern.selective_counter, memval)
    # apply the memory onto the POMDP to get the cartesian product
    pomdp = stormpy.pomdp.unfold_memory(pomdp, memory)
    print("Number of pomdp states before simple:",pomdp.nr_states)
    print("Number of transitions: {}".format(pomdp.nr_transitions))

    # make the POMDP simple. This step is optional but often beneficial
    pomdp = stormpy.pomdp.make_simple(pomdp)
    print("Number of pomdp states after simple:",pomdp.nr_states)
    print("Number of transitions: {}".format(pomdp.nr_transitions))

    # apply the unknown FSC to obtain a pmc from the POMDP
    pmc = stormpy.pomdp.apply_unknown_fsc(pomdp, stormpy.pomdp.PomdpFscApplicationMode.simple_linear)
    print("Number of pomdp states after simple",pmc.nr_states)
    print("Number of transitions: {}".format(pmc.nr_transitions))
    print("applied pmc")
    path_pmc = "export_" + str(memval) + "_mem_" + model_path
    stormpy.export_to_drn(pmc, path_pmc)
    print("built model")
    #print(model.initial_states)
    fsc_parameters = pmc.collect_probability_parameters() - pomdp.collect_probability_parameters()
    print("number of pomdp parameters:",len(fsc_parameters))
    #print(fsc_parameters)
    #print(pomdp_parameters)
    # gather intervals and items



    print("Processing intervals\n")

    intervals, polyhedrons, items = interval_parser.parse_model_interval(pmc,interval_parameters,interval_path)

   # print("\n\n")
    print("Computing vertices")

    for p in polyhedrons:
        #print(p.vertices)
        p.compute_vertices()
        #print(p.vertices)
    #    print(type(p.state))
    t1 = time.time()
    print("Vertex computation time: ", str(t1 - t0))


    #for item in items:
    #    print(item,"printing in the main")

    properties = stormpy.parse_properties(formula_str)
    print("Building model from {}".format(model_path))
    #parameters = model.collect_probability_parameters()
    print(pmc.nr_states)
    #compute prob01max states
    prob0E, prob1A = stormpy.prob01max_states(pmc, properties[0].raw_formula.subformula)

    #threshold = 0.90
    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=0.05, maxiter=1000, graph_epsilon=1e-2, silent=False)

    # run solver

    solver = QcqpSolver()

   # for x in fsc_parameters:
   #     print(x)
   #     print(type(x))
   #     print(x.id)
   #     print(type(x.id))
   #     break

    fsc_parameters_ids = set()
    for x in fsc_parameters:
        fsc_parameters_ids.add(x.id)
    interval_parameters_ids = set()
    for x in interval_parameters:
        interval_parameters_ids.add(x.id)

    print("preparing polyhedrons")
    prep0 = time.time()
    polyhedron_state_map = dict()
    for id in range(pmc.nr_states):
        current_polyhedrons = []
        for p in polyhedrons:
            if int(p.state) == int(id):
                current_polyhedrons.append(p)
        polyhedron_state_map[id] = current_polyhedrons
    prep1 = time.time()
    print("Polyhedron time: " + str(prep1 - prep0))
    print("Total preparation time: ", str(prep1 - t0))
    print("Starting solver")
    result = solver.run(pmc, fsc_parameters, fsc_parameters_ids, interval_parameters, interval_parameters_ids, properties, prob0E, prob1A, threshold, direction, options, intervals, polyhedrons, polyhedron_state_map, True)

    t_end = time.time()

    print("number of iterations={}".format(solver.iterations))
    print("Vertex computation time: ", str(t1 - t0))
    print("Polyhedron time: " + str(prep1 - prep0))
    print("Total preparation time: ", str(prep1 - t0))
    print("solver time={}".format(solver.solver_timer))
    print("total time={}".format(str(t_end - t0)))

    #compute the policy against the robust interval
#    interval_path="satellite_robust.intervals"
#    intervals, polyhedrons,items = interval_parser.parse_model_interval(pmc,pomdp_parameters,interval_path)
#    regiondict = dict()
#    for x in pomdp_parameters:
#        for item in items:
#            if item.name == x.name:
#                regiondict[x] = (stormpy.RationalRF(item.lowerbound), stormpy.RationalRF(item.upperbound))
#    region = stormpy.pars.ParameterRegion(regiondict)
#    instantiator = stormpy.pars.PartialPDtmcInstantiator(pmc)
#    instantiated_model = instantiator.instantiate(solver.solver_params)
                    #assert interval_parameters == instantiated_model.collect_probability_parameters()

#    env = stormpy.Environment()
#    env.solver_environment.set_linear_equation_solver_type(stormpy.EquationSolverType.eigen)
#    env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.policy_iteration
                    #env.solver_environment.minmax_solver_environment.method = stormpy.MinMaxMethod.optimistic_value_iteration

#    start_check = time.time()

#    region_checker = stormpy.pars.create_region_checker(env, instantiated_model, properties[0].raw_formula,
#                                                                        allow_model_simplification=False)

#    print("region check")
#    result = region_checker.get_bound_all_states(env, region, maximise=False)
#    end_check = time.time()
#    print("model checking ans:")
#    ansval=result.at(pmc.initial_states[0])
#    print(ansval)





# small debug thing to figure out how things are stored
def test():
    model_path = "4x4_grid_trap_interval_multitrap.drn"
    interval_path = "4x4_grid_trap_long.intervals"

    formula_str = 'P=? [F "goal" ]'

    # adapt into lower and upper model
    model = stormpy.build_parametric_model_from_drn(model_path)

    print(model.initial_states)

    intervals, polyhedrons = interval_parser.parse_input(interval_path)
    for p in polyhedrons:
        p.compute_vertices()


    properties = stormpy.parse_properties(formula_str)
    print("Building model from {}".format(model_path))

    threshold = 0.5
    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=0.05, maxiter=10000, graph_epsilon=1e-6, silent=False)

    parameters = model.collect_probability_parameters()
    prob0E, prob1A = stormpy.prob01max_states(model, properties[0].raw_formula.subformula)
    print(prob1A)
    #

    # result = solver.run(reward_model_name ,model_rew, parameters_rew, rew0, rew_threshold, direction, options)
    solver = QcqpSolver()
    result = solver.run(model, parameters, properties, prob0E, prob1A, threshold, direction, options, polyhedrons)
    #print("RESULT: " + str(result))

    #parameter_assignments = dict([[x, result.parameter_values[x.id]] for x in parameters if x.name[0] != 'I'])
    #print(parameter_assignments)

    #instantiator = stormpy.pars.PDtmcInstantiator(model)

    # # Check distance to result by storm (notice that also the point that storm checks slightly differs)
    #rational_parameter_assignments = dict([[x, stormpy.RationalRF(val)] for x, val in parameter_assignments.items()])
    # print(rational_parameter_assignments)
    #instantiated_model = instantiator.instantiate(rational_parameter_assignments)
    #mc_res = stormpy.model_checking(instantiated_model, properties[0]).at(model.initial_states[0])

    #print("Qcqp: {}: Mc: {}".format(result.value_at_initial, mc_res))
    print("iternum={}".format(solver.iterations))
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


    threshold = 0.7
    direction = "below"  # can be "below" or "above"
    options = QcqpOptions(mu=0.05, maxiter=10000, graph_epsilon=1e-6, silent=False)

    parameters = model.collect_probability_parameters()
    prob0E, prob1A = stormpy.prob01max_states(model, properties[0].raw_formula.subformula)
    print(prob1A)
    #

    #result = solver.run(reward_model_name ,model_rew, parameters_rew, rew0, rew_threshold, direction, options)
    solver = QcqpSolver()
    result = solver.run(model, parameters, properties, prob0E, prob1A, threshold, direction, options, polyhedrons)
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
    model_path, interval_path, formula_str, threshold = input_files()
    main(model_path, interval_path, formula_str, threshold)


