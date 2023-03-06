import os
import stormpy
import stormpy.pomdp
import stormpy.pars
from pycarl.cln.cln import Rational
import time

class PMC:

    def __init__(self, model):
        self.model = model
        self.nr_parameters = len(self.model.collect_probability_parameters())
        self.result = None  
    
    def instantiate_parameters(self, parameter_values):
        parameters = list(self.model.collect_probability_parameters())
        if len(parameter_values) == len(parameters):
            instantiator = stormpy.pars.PDtmcInstantiator(self.model)
            par_dict = dict()
            for x in range(len(parameters)):
                parameter_x = parameters[x]
                par_dict[parameter_x] = Rational(parameter_values[x])
            self.model = instantiator.instantiate(par_dict)
        else:
            raise Exception("Number of parameter values does not match number of parameters, expected:"\
                , len(parameters), ", found:", len(parameter_values)) 

    # Perform model checking using a property
    def model_checking(self, property):
        self.inspect_model()
        self.result = stormpy.model_checking(self.model, property)
    
    def get_results(self):
        #Qualitative result types do not have a get_values() method 
        #but instead a get_truth_values() method
        qualitative_result_types = [stormpy.ExplicitQualitativeCheckResult, \
            stormpy.SymbolicQualitativeCheckResult]
        
        #These result type do not have a get_values() or a get_truth_values() method
        not_supported_result_types = [stormpy.SymbolicParametricQuantitativeCheckResult, \
            stormpy.SymbolicExactQuantitativeCheckResult, stormpy.SymbolicQuantitativeCheckResult]

        check_result_type = lambda type: isinstance(self.result, type)

        if any(map(check_result_type, not_supported_result_types)):
            raise Exception("Result type not supported")
        elif any(map(check_result_type, qualitative_result_types)):
            return self.result.get_truth_values()
        else:
            return self.result.get_values() 

    # Inspect model for invalid or missing components
    def inspect_model(self):
        if not self.has_transition_matrix():
            raise Exception("The pMC model does not have a transition matrix")

        if not self.has_state_labels():
            print("No state labels where found in the pMC model")

        if not self.has_reward_model():
            raise Exception("The pMC model is missing reward models")

        if not self.model.has_state_valuations():
            print("The pMC model has no state valuations")
        
        if self.has_parameters():
            raise Exception("The pMC model has uninstantiated parameters")
    
    def has_transition_matrix(self):
        return self.model.nr_states > 0 \
            and self.model.nr_transitions > 0 \
            and len(self.model.transition_matrix) > 0

    def has_state_labels(self):
        return len(self.model.labeling.get_labels()) > 0

    def has_reward_model(self):
        return True if len(self.model.reward_models) > 0 else False

    def has_parameters(self):
        return self.model.has_parameters

    # Get functions
    def get_model_type(self):
        return self.model.model_type
    
    def get_states(self):
        return self.model.states
    
    def get_nr_states(self):
        return len(self.model.states)
    
    def get_state_labels(self):
        return self.model.labeling.get_labels()
    
    def get_reward_models(self):
        return self.model.reward_models

    # Print functions
    def print_model_type(self):
        print("Model type:",self.model.model_type)

    def print_nr_states(self):
        print("Number of states:", self.model.nr_states)
    
    def print_actions(self):
        for state in self.model.states:
            print("State", state, " has actions", [str(action) for action in state.actions])

    def print_transitions(self):
        for state in self.model.states:
            for action in state.actions:
                for transition in action.transitions:
                    print("State", state, " using action", action," goes to state",\
                            transition.column," with propability", transition.value())
    
    def print_rewards(self):
        for key in self.model.reward_models:
            reward_model = self.model.reward_models[key]
            if reward_model.has_state_rewards:
                print("has state rewards")
                self.print_state_reward(reward_model.state_rewards)
            if reward_model.has_state_action_rewards:
                self.print_state_action_reward(reward_model.state_action_rewards)
            if reward_model.has_transition_rewards:
                print("has transition rewards")
    
    def print_state_reward(self, reward):
        print("***** State Rewards *****")
        i = 0
        for state in self.get_states():
            for action in state.actions:
                print("State", state, " has a reward of", reward[i])
                i += 1


    def print_state_action_reward(self, reward):
        print("***** State Action Rewards *****")
        i = 0
        for state in self.get_states():
            for action in state.actions:
                print("State", state, " using action", action, " has a reward of", reward[i])
                i += 1

    def print_results(self):
        cnt = 0
        for x in self.get_results():
            print("Result", cnt," =", x)
            cnt += 1
    
    def export_parametric_to_drn(self, file):
        stormpy.export_parametric_to_drn(self.model, file)
