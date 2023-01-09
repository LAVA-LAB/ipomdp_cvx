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
import interval_parser

def value_iteration_iterator(model,parameters, solution, interval_parameters_ids ,prob0E,prob1A,polyhedrons,option):
    point = dict()

    for p in parameters:
        point[p.id] = solution[p.id]

    numstate = model.nr_states

    induced_mc_nominal = [[0 for _ in range(numstate)] for __ in range(numstate)]
    storedprobs=[1  for _ in range(numstate)]
    storedprobsold=[1  for _ in range(numstate)]
    numiter=50
    for iter in range(numiter):

        for state in model.states:
            if prob1A.get(state):
                storedprobs[state.id]=1
                storedprobsold[state.id]=1
                #print("lel")
            if prob0E.get(state):
                storedprobs[state.id]=0
                storedprobsold[state.id]=0

        for state in model.states:
            # print(state.id)
            # print (prob1.get(state))
            # Cons=values constraints on the right hand side for a pdtmc

            # Find the polyhedrons at the current state

            current_polyhedrons = []
            for p in polyhedrons:
                if int(p.state) == int(state.id):
                    current_polyhedrons.append(p)


            if len(current_polyhedrons) == 0:
                probval=0

                for action in state.actions:
                    for transition in action.transitions:
                        transition_value = transition.value()
                        succ = int(transition.column)
                        # Value of transition
                        transition_value = transition.value()
                        # TODO check that the denominator is indeed constant?
                        # Denominator of transition
                        den = transition_value.denominator.constant_part()
                        denom1 = 1 / float(den)
                        if not transition_value.is_constant():
                            num = transition_value.numerator.polynomial()
                            for t in num:
                                # If the transition term is a constant
                                if t.is_constant():
                                    # Add just value of transition is the successor state is prob1 to the constraints
                                    if prob1A.get(succ):
                                        probval+= float(t.coeff) * denom1
                                    # Add nothing successor state is prob0
                                    elif prob0E.get(succ):
                                        pass
                                    # Else add transitionvalue*p_succ to the constraint
                                    else:
                                        probval+= float(t.coeff) * denom1 * storedprobsold[succ]
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    param_id = t.monomial[0][0].id
                                    coeff = float(t.coeff)
                                    # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                    if prob1A.get(succ):
                                        probval+= coeff * float(point[param_id]) * denom1
                                    # Add nothing successor state is prob0
                                    elif prob0E.get(succ):
                                        pass
                                    # Add the quadratic term to the constraints
                                    else:
                                        probval+= coeff * float(point[param_id]) * denom1* storedprobsold[succ]
                storedprobs[state.id]=probval

            else:
                # robust constraint
                # get all possible instantiations
                combinations = interval_parser.calculate_vertex_combinations(current_polyhedrons)
                # for each instantiation build a constraint
                valueiterval=1e8
                for instantiation in combinations:
                    # for every instantiation, a constraint is added
                    probval = 0
                    #print(state.id)


                    for action in state.actions:
                        for transition in action.transitions:

                            transition_value = transition.value()
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
                                            probval+= float(t.coeff) * denom1
                                        # Add nothing successor state is prob0
                                        elif prob0E.get(succ):
                                            pass
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            probval+= float(t.coeff) * denom1 * storedprobsold[succ]

                                    # If the transition term is not a constant
                                    else:
                                        # if t.tdeg > 1:
                                        #    raise RuntimeError("We expect the term to be a single variable")
                                        # print(t)
                                        if t.tdeg == 1:
                                            # term of degree 1, check if it is an interval or a real parameter
                                            param_name = t.monomial[0][0].name
                                            param_id = t.monomial[0][0].id
                                            if param_id in interval_parameters_ids:
                                                #if param_name[0] == 'I':

                                                interval_name = param_name
                                                interval_value = -1
                                                for vertex in instantiation:
                                                    for vertex_value in vertex.vertex_values:
                                                        if vertex_value.name == interval_name:
                                                            interval_value = vertex_value.value
                                                            coeff = float(t.coeff) * float(interval_value)
                                                            break
                                                if interval_value == -1:
                                                    raise RuntimeError("Interval not part of the instantiation")
                                                probval+= coeff * denom1 * storedprobsold[succ]
                                                param_id = None

                                            else:
                                                param_id = t.monomial[0][0].id
                                                coeff = float(t.coeff)
                                                probval+= coeff * denom1 * storedprobsold[succ] * float(point[param_id])


                                        elif t.tdeg == 2:
                                            par1 = t.monomial[0][0].id
                                            par2 = t.monomial[1][0].id
                                            if par1 in interval_parameters_ids:
                                                # print(par1)
                                                param_id = t.monomial[1][0].id
                                                interval_id = par1
                                            elif par2 in interval_parameters_ids:
                                                # print(par2)
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
                                                        coeff = float(t.coeff) * float(interval_value)
                                                        break
                                            if interval_value == -1:
                                                print("interval id: " + str(interval_id))
                                                raise RuntimeError("Interval not part of the instantiation")
                                            #print(interval_value,coeff,interval_name,point[param_id])
                                            probval += coeff * denom1 * storedprobsold[succ] * float(
                                                point[param_id])


                                        else:
                                            raise RuntimeError("We expect the term to be a single variable")
         #           print(probval,valueiterval,state.id)
                    if option=="min":
                        if valueiterval>probval:
                            storedprobs[state.id]=probval
                            valueiterval=probval
                    elif option=="max":
                        if valueiterval<probval:
                            storedprobs[state.id]=probval
                            valueiterval=probval
                    else:
                        raise RuntimeError("Expecting either min or max as option")
        for state in model.states:
            storedprobsold[state.id]=storedprobs[state.id]

        print("Final reachability probs:")
        print(storedprobs)

    # only return reachability at initial state
    return storedprobs[model.initial_states[0]]

def value_iteration_rew_iterator(model, parameters, rew0, polyhedrons,option):
    reward_model_name = list(model.reward_models.keys())[0]
    point = dict()

    for p in parameters:
        param_string = str(p)
        if not param_string[0] == 'I':
            point[p.id] = 0.3

    numstate = model.nr_states

    induced_mc_nominal = [[0 for _ in range(numstate)] for __ in range(numstate)]
    storedprobs = [0 for _ in range(numstate)]
    storedprobsold = [0 for _ in range(numstate)]
    numiter = 50
    for iter in range(numiter):

        for state in model.states:
            if rew0.get(state):
                storedprobs[state.id] = 0
                storedprobsold[state.id] = 0
                # print("lel")


        for state in model.states:
            # print(state.id)
            # print (prob1.get(state))
            # Cons=values constraints on the right hand side for a pdtmc

            # Find the polyhedrons at the current state

            current_polyhedrons = []
            for p in polyhedrons:
                if p.state == str(state.id):
                    current_polyhedrons.append(p)
            reward_at_state = model.reward_models[reward_model_name].state_rewards[int(state.id)]
            reward_at_state=float(str(reward_at_state))

            if len(current_polyhedrons) == 0:
                probval = 0


                for action in state.actions:
                    for transition in action.transitions:
                        transition_value = transition.value()
                        succ = int(transition.column)
                        # Value of transition
                        transition_value = transition.value()
                        # TODO check that the denominator is indeed constant?
                        # Denominator of transition
                        den = transition_value.denominator.constant_part()
                        denom1 = 1 / float(den)
                        if not transition_value.is_constant():
                            num = transition_value.numerator.polynomial()
                            for t in num:
                                # If the transition term is a constant
                                if t.is_constant():
                                    # Add just value of transition is the successor state is prob1 to the constraints
                                    #if prob1A.get(succ):
                                        #probval += float(t.coeff) * denom1
                                    # Add nothing successor state is prob0
                                    if rew0E.get(succ):
                                        pass
                                    # Else add transitionvalue*p_succ to the constraint
                                    else:
                                        probval += reward_at_state + float(t.coeff) * denom1 * storedprobsold[succ]
                                else:
                                    if t.tdeg > 1:
                                        raise RuntimeError("We expect the term to be a single variable")
                                    param_id = t.monomial[0][0].id
                                    coeff = float(t.coeff)
                                    # Adds transitionvalue*parameter_variable to the constraint if the successor is prob1

                                    #if prob1A.get(succ):
                                        #probval += coeff * float(point[param_id]) * denom1
                                    # Add nothing successor state is prob0
                                    if rew0E.get(succ):
                                        pass
                                    # Add the quadratic term to the constraints
                                    else:
                                        probval += reward_at_state + coeff * float(point[param_id]) * denom1 * storedprobsold[succ]
                storedprobs[state.id] = probval

            else:
                # robust constraint
                # get all possible instantiations
                combinations = interval_parser.calculate_vertex_combinations(current_polyhedrons)
                # for each instantiation build a constraint
                valueiterval = 1e8
                for instantiation in combinations:
                    # for every instantiation, a constraint is added
                    probval = 0
                    # print(state.id)

                    for action in state.actions:
                        for transition in action.transitions:

                            transition_value = transition.value()
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
                                        #if prob1A.get(succ):
                                            #probval += float(t.coeff) * denom1
                                        # Add nothing successor state is prob0
                                        if rew0.get(succ):
                                            pass
                                        # Else add transitionvalue*p_succ to the constraint
                                        else:
                                            probval += reward_at_state + float(t.coeff) * denom1 * storedprobsold[succ]

                                    # If the transition term is not a constant
                                    else:
                                        # if t.tdeg > 1:
                                        #    raise RuntimeError("We expect the term to be a single variable")
                                        # print(t)
                                        if t.tdeg == 1:
                                            # term of degree 1, check if it is an interval or a real parameter
                                            param_name = t.monomial[0][0].name
                                            if param_name[0] == 'I':
                                                param_id = None
                                                interval_name = param_name
                                                interval_value = -1
                                                for vertex in instantiation:
                                                    for vertex_value in vertex.vertex_values:
                                                        if vertex_value.name == interval_name:
                                                            interval_value = vertex_value.value
                                                            coeff = float(t.coeff) * float(interval_value)
                                                            break
                                                if interval_value == -1:
                                                    raise RuntimeError("Interval not part of the instantiation")
                                                probval += coeff * denom1 * storedprobsold[succ]

                                            else:
                                                param_id = t.monomial[0][0].id
                                                coeff = float(t.coeff)
                                                probval += coeff * denom1 * storedprobsold[succ] * float(
                                                    point[param_id])


                                        elif t.tdeg == 2:
                                            par1 = t.monomial[0][0].name
                                            par2 = t.monomial[1][0].name
                                            if par1[0] == 'I':
                                                # print(par1)
                                                param_id = t.monomial[1][0].id
                                                interval_name = par1
                                            elif par2[0] == 'I':
                                                # print(par2)
                                                param_id = t.monomial[0][0].id
                                                interval_name = par2
                                            else:
                                                raise RuntimeError(
                                                    "Multiple non-interval parameters in the same monomial")

                                            # find value of the interval given instantiation
                                            interval_value = -1
                                            for vertex in instantiation:
                                                for vertex_value in vertex.vertex_values:
                                                    if vertex_value.name == interval_name:
                                                        interval_value = vertex_value.value
                                                        coeff = float(t.coeff) * float(interval_value)
                                                        break
                                            if interval_value == -1:
                                                print("interval name: " + str(interval_name))
                                                raise RuntimeError("Interval not part of the instantiation")
                                            # print(interval_value,coeff,interval_name,point[param_id])
                                            probval += reward_at_state + coeff * denom1 * storedprobsold[succ] * float(
                                                point[param_id])


                                        else:
                                            raise RuntimeError("We expect the term to be a single variable")
                    #           print(probval,valueiterval,state.id)
                    if option=="min":
                        if valueiterval>probval:
                            storedprobs[state.id]=probval
                            valueiterval=probval
                    elif option=="max":
                        if valueiterval<probval:
                            storedprobs[state.id]=probval
                            valueiterval=probval
                    else:
                        raise RuntimeError("Expecting either min or max as option")

        for state in model.states:
            storedprobsold[state.id] = storedprobs[state.id]

        print("Final reachability probs:")
        print(storedprobs)

    # only return reachability at initial state
    return storedprobs[model.initial_states[0]]

    #pass


# small debug thing to figure out how things are stored
def test():
    model_path = "input/pMC1.drn"
    interval_path = "input/pMC1.intervals"

    formula_str = 'P=? [F "goal" ]'

    # adapt into lower and upper model
    model = stormpy.build_parametric_model_from_drn(model_path)

    #print(model.initial_states)

    intervals, polyhedrons = interval_parser.parse_input(interval_path)
    for p in polyhedrons:
        p.compute_vertices()


    properties = stormpy.parse_properties(formula_str)
    print("Building model from {}".format(model_path))
    parameters = model.collect_probability_parameters()
    prob0E, prob1A = stormpy.prob01max_states(model, properties[0].raw_formula.subformula)
    value_iteration_iterator(model,parameters,prob0E,prob1A,polyhedrons,"min")

def test_rew():
    model_path = "input/pMC1_rew.drn"
    interval_path = "input/pMC1.intervals"

    formula_str = 'R=? [F "goal" ]'

    # adapt into lower and upper model
    model = stormpy.build_parametric_model_from_drn(model_path)

    #print(model.initial_states)

    intervals, polyhedrons = interval_parser.parse_input(interval_path)
    for p in polyhedrons:
        p.compute_vertices()


    properties = stormpy.parse_properties(formula_str)
    print("Building model from {}".format(model_path))
    parameters = model.collect_probability_parameters()
    rew0 = find_rew0_states(model, properties[0])

    value_iteration_rew_iterator(model,parameters,rew0,polyhedrons,"min")


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


def find_rew0_states(model,property):
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
    test_rew()


