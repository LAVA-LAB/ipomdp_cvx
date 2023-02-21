## @package pomdp_solver
#  This is a script developed to automate the calculation of a policy satisfying a property 
#  for a given pomdp model. 
#
#  This script contains the solve_pomdp function for an automated way to calculate a policy. The script
#  also contains the POMDP and PMC classes that provide functionalities to explore the POMDP and 
#  pMC models respectively. 
#
#  This script requires the following tools: Storm, Stormpy, Carl, Pycarl, Gurobi
#
#  This script uses the solve_nominal_pomdp script provided by Nils Jansen (n.jansen@science.ru.nl) 
#  and Marnix Suilen (marnix.suilen@ru.nl) of the Radboud Univeristy.
#
#  Created by Anass Fakir (a.fakir@student.ru.nl).

import os
import stormpy
import stormpy.pomdp
import stormpy.pars
from pycarl.cln.cln import Rational
import time
import solver
import pycarl

## Function to automate the calculation of a policy for a pomdp model using Storm and Gurobi. 
#
#  @param model Path to a pomdp model
#  @param properties String of a property that the calculated policy must satisfy.
#  @param memory_size Integer for the memory size that is unfolded in the pomdp model.
#  @param drn_export_file File to which the pMC model will be written in the DRN format.
#  @param export_instantiated_model Boolean to indicate whether to export the instantiated pMC model.
def solve_pomdp(model, properties, memory_size=2, 
                drn_export_file="temp.drn", 
                export_instatiated_model=False):
    pomdp = POMDP()
    pomdp.load_model(model)
    pomdp.build_model()

    ##Unfold memory and apply a fsc to get an pMC
    pmc = PMC(pomdp.build_pmc(nr_memory_state=memory_size))
    ##Export the pMC to a DRN file
    pmc.export_to_drn(drn_export_file)
    ##Use the solver script to obtain parameters for the pMC
    result = solver.solve_nominal_pomdp(drn_export_file, properties, 0.99, silent=False)
    if export_instatiated_model:
        ##Instantiate the pMC model
        pmc.instantiate_parameters(list(result.parameter_values.values()))
        ##Export the instantiated pMC model.
        pmc.export_to_drn("instantiated_" + drn_export_file)
    return result
    

## Class to build and explore a POMDP model.
#
#  Supports DRN and Prism models. The use of Jani files has not been tested to give any guarantees. 
class POMDP:
    ## The constuctor
    def __init__(self):
        ## Can be drn, prism or jani.
        self.file_type = ""
        ## Path to the model.
        self.file_path = ""
        self.properties = None
        ## Parsed pomdp model
        self.program = None
        ## Build pomdp model
        self.model = None
        
    ## Load a POMDP model.
    #
    #  @param model_file Path to a pomdp model
    def load_model(self, model_file):
        self.properties = None
        ## Obtain the file type from the path.
        self.file_type = os.path.splitext(model_file)[1]
        self.file_path = model_file
        ## Parse the pomdp model.
        if self.file_type == ".prism":
            self.program = stormpy.parse_prism_program(model_file)
        elif self.file_type == ".drn":
            self.program = None
        elif self.file_type == ".jani":
            self.program, self.properties = stormpy.parse_jani_model(model_file)
        else:
            raise Exception("Model file type not supported") 
    
    ## Parse property string.
    #
    #  @param properties Property string
    def parse_properties(self, properties):
        if self.file_type == ".prism":
            self.properties = stormpy.parse_properties_for_prism_program(properties, self.program)
        elif self.file_type == ".jani":
            self.properties = stormpy.parse_properties_for_jani_model(properties, self.program)
        elif self.file_type == ".drn":
            self.properties = stormpy.parse_properties_without_context(properties)
        else:
            raise Exception("Model file type not supported to parse properties")      

    ## Instantiate undefined constants of a pomdp model.
    #
    #  @param constants a String containing a variable name and its value(e.g. "K=4, T=2")
    def set_undefined_constants(self, constants):
        desc = stormpy.SymbolicModelDescription(self.program)
        #Parsing the constants
        const = desc.parse_constant_definitions(constants)
        #Instantiate the constants
        if self.file_type == ".prism":
            self.program = desc.instantiate_constants(const).as_prism_program()
        elif self.file_type == ".jani":
            self.program = desc.instantiate_constants(const).as_jani_model()
        else:
            raise Exception("Model file type not supported to define constants")    

    ## Build the POMDP model.
    #
    def build_model(self):
        ## Different build functions depending on the model format.
        if self.file_type == ".prism" or self.file_type == ".jani":
            ## Defining build options to ensure that choice labels remain.
            options = stormpy.BuilderOptions()
            options.set_build_state_valuations()
            options.set_build_choice_labels()
            self.model = stormpy.build_sparse_model_with_options(self.program, options)
        elif self.file_type == ".drn":
            ## Defining build options to ensure that choice labels remain.
            parser_options = stormpy.DirectEncodingParserOptions()
            parser_options.build_choice_labels = True
            self.model = stormpy.build_model_from_drn(self.file_path, parser_options)   
        else:
            raise Exception("Model file type not supported to build models") 

    ## Unfold memory into the POMDP model and apply a FSC to obtain a pMC.
    #
    #  @param nr_memory_state Number of memory nodes.
    #  @param pomdp_memory_pattern Memory pattern.
    #  @param fsc_application_mode Application mode for the FSC.
    #  @return A pMC model.
    def build_pmc(self,\
        nr_memory_state,\
        pomdp_memory_pattern=stormpy.pomdp.PomdpMemoryPattern.selective_counter,\
        fsc_application_mode=stormpy.pomdp.PomdpFscApplicationMode.simple_linear):
        ## Make the pomdp model canonical.
        self.model = stormpy.pomdp.make_canonic(self.model)

        memoryBuilder = stormpy.pomdp.PomdpMemoryBuilder()
        memory = memoryBuilder.build(pomdp_memory_pattern, nr_memory_state)
        ## Unfold memory into the POMDP model.
        self.model = stormpy.pomdp.unfold_memory(self.model, memory, add_memory_labels=True, \
        keep_state_valuations=True)
        ## Make POMDP model simple.
        self.model = stormpy.pomdp.make_simple(self.model, keep_state_valuations=True)
        ## Apply a FSC to obtain a pMC model.
        pmc = stormpy.pomdp.apply_unknown_fsc(self.model, fsc_application_mode)
        return pmc

    ## Inspect the POMDP model for invalid or missing components.
    #
    def inspect_model(self):
        if not self.has_transition_matrix():
            raise Exception("The POMDP model does not have a transition matrix")

        if not self.has_state_labels():
            print("No state labels where found in the POMDP model")

        if not self.has_reward_model():
            raise Exception("The PODP model is missing reward models")
        '''
        if not self.model.has_choice_labeling():
            raise Exception("The POMDP model does not have choice labels")
        '''
        if not self.model.has_state_valuations():
            print("The POMDP model has no state valuations")
    
    ## Check if the POMDP model has a transition matrix.
    #
    #  @return Boolean.
    def has_transition_matrix(self):
        return self.model.nr_states > 0 \
            and self.model.nr_transitions > 0 \
            and len(self.model.transition_matrix) > 0

    ## Check if the POMDP model has a state labels.
    #
    #  @return Boolean.
    def has_state_labels(self):
        return len(self.model.labeling.get_labels()) > 0

    ## Check if the POMDP model has a reward model.
    #
    #  @return Boolean.
    def has_reward_model(self):
        return True if len(self.model.reward_models) > 0 else False
    
    ## Check if the POMDP model has undefined constants.
    #
    #  @return Boolean.
    def has_undefined_constants(self):
        return self.program.has_undefined_constants

    ## Get the POMDP model type.
    #
    #  @return Model type.
    def get_model_type(self):
        return self.model.model_type

    ## Get the states of the POMDP model.
    #
    #  @return List of states.
    def get_states(self):
        return self.model.states
    
    ## Get the number of states in the POMDP model.
    #
    #  @return Number of states.
    def get_nr_states(self):
        return len(self.model.states)

    ## Get the state labels of the POMDP model.
    #
    #  @return List of state labels.
    def get_state_labels(self):
        return self.model.labeling.get_labels()
    
    '''
    def get_choice_labels(self):
        return self.model.choice_labeling.get_labels()
    '''

    ## Get the reward model of the POMDP model.
    #
    #  @return reward model.
    def get_reward_models(self):
        return self.model.reward_models
    
    ## Get the observations of the POMDP model.
    #
    #  @return List of observations.
    def get_observations(self):
        return self.model.observations

    ## Print the model type.
    #
    def print_model_type(self):
        print("Model type:", self.get_model_type())

    ## Print the number of states.
    #
    def print_nr_states(self):
        print("Number of states:", self.get_nr_states())
    
    ## Print the actions.
    #
    def print_actions(self):
        for state in self.model.states:
            print("State", state, " has actions", [str(action) for action in state.actions])

    ## Print the transitions of the model.
    #
    def print_transitions(self):
        for state in self.model.states:
            for action in state.actions:
                for transition in action.transitions:
                    print("State", state, " using action", action," goes to state",\
                         transition.column," with propability", transition.value())
    
    ## Print the observations of the model.
    #
    def print_observations(self):
        for i in range(self.get_nr_states()):
            print("State", i, " has observation", self.model.get_observation(i))
    
    ## Print the reward model.
    #
    def print_rewards(self):
        for key in self.model.reward_models:
            reward_model = self.model.reward_models[key]
            if reward_model.has_state_rewards:
                self.print_state_reward(reward_model.state_rewards)
            if reward_model.has_state_action_rewards:
                self.print_state_action_reward(reward_model.state_action_rewards)
            if reward_model.has_transition_rewards:
                print("has transition rewards")

    ## Print the state rewards.
    #
    def print_state_reward(self, reward):
        print("***** State Rewards *****")
        i = 0
        for state in self.get_states():
            for action in state.actions:
                print("State", state, " has a reward of", reward[i])
                i += 1

    ## Print the state action rewards.
    #
    def print_state_action_reward(self, reward):
        print("***** State Action Rewards *****")
        i = 0
        for state in self.get_states():
            for action in state.actions:
                print("State", state, " using action", action, " has a reward of", reward[i])
                i += 1
    
    ## Print constants of the model.
    #
    def print_constants(self):
        for c in self.program.constants:
            if c.defined:
                print("Constant:", c.name, " type:", c.type, " Definition:", c.definition)
            else:
                print("Constant:", c.name, " type:", c.type, " Definition: Undefined")

## Class to explire and export a pMC model.
#
#  Takes a pMC model as an argument that is returned by the buil_pmc function of the POMDP class. 
class PMC:
    ## The constuctor
    def __init__(self, model):
        ## The pMC model.
        self.model = model
        ## Number of parameters of the pMC model.
        self.nr_parameters = len(self.model.collect_probability_parameters())
    
    ## Instantiate a pMC.
    #
    #  @param parameter_values List of values to instantiate parameters.
    def instantiate_parameters(self, parameter_values):
        ## Collect the parameters.
        parameters = list(self.model.collect_probability_parameters())
        ## Check if the number of parameters equals the number of values.
        if len(parameter_values) == len(parameters):
            instantiator = stormpy.pars.PDtmcInstantiator(self.model)
            ## Create a dictionary of parameters to values
            par_dict = dict()
            for x in range(len(parameters)):
                parameter_x = parameters[x]
                par_dict[parameter_x] = Rational(parameter_values[x])
                print("Instantiated parameter ", str(parameter_x), "with value ", str(parameter_values[x]))
            ## Instatiate the model
            self.model = instantiator.instantiate(par_dict)
        else:
            raise Exception("Number of parameter values does not match number of parameters, expected:"\
                , len(parameters), ", found:", len(parameter_values)) 

    ## Inspect model for invalid or missing components
    #
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
    
    ## Check if the pMC model has a transition matrix.
    #
    #  @return Boolean.
    def has_transition_matrix(self):
        return self.model.nr_states > 0 \
            and self.model.nr_transitions > 0 \
            and len(self.model.transition_matrix) > 0

    ## Check if the pMC model has a state labels.
    #
    #  @return Boolean.
    def has_state_labels(self):
        return len(self.model.labeling.get_labels()) > 0

    ## Check if the pMC model has a reward model.
    #
    #  @return Boolean.
    def has_reward_model(self):
        return True if len(self.model.reward_models) > 0 else False

    ## Check if the pMC model has uninstantiated parameters.
    #
    #  @return Boolean.
    def has_parameters(self):
        return self.model.has_parameters

    ## Get the pMC model type.
    #
    #  @return Model type.
    def get_model_type(self):
        return self.model.model_type
    
    ## Get the states of the pMC model.
    #
    #  @return List of states.
    def get_states(self):
        return self.model.states
    
    ## Get the number of states in the pMC model.
    #
    #  @return Number of states.
    def get_nr_states(self):
        return len(self.model.states)
    
    ## Get the state labels of the pMC model.
    #
    #  @return List of state labels.
    def get_state_labels(self):
        return self.model.labeling.get_labels()
    
    ## Get the reward model of the pMC model.
    #
    #  @return reward model.
    def get_reward_models(self):
        return self.model.reward_models

    ## Print the model type.
    #
    def print_model_type(self):
        print("Model type:",self.model.model_type)

    ## Print the number of states.
    #
    def print_nr_states(self):
        print("Number of states:", self.model.nr_states)
    
    ## Print the actions.
    #
    def print_actions(self):
        for state in self.model.states:
            print("State", state, " has actions", [str(action) for action in state.actions])

    ## Print the transitions of the model.
    #
    def print_transitions(self):
        for state in self.model.states:
            for action in state.actions:
                for transition in action.transitions:
                    print("State", state, " using action", action," goes to state",\
                            transition.column," with propability", transition.value())
    
    ## Print the reward model.
    #
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
    
    ## Print the state rewards.
    #
    def print_state_reward(self, reward):
        print("***** State Rewards *****")
        i = 0
        for state in self.get_states():
            for action in state.actions:
                print("State", state, " has a reward of", reward[i])
                i += 1

    ## Print the state action rewards.
    #
    def print_state_action_reward(self, reward):
        print("***** State Action Rewards *****")
        i = 0
        for state in self.get_states():
            for action in state.actions:
                print("State", state, " using action", action, " has a reward of", reward[i])
                i += 1

    ## Export the pMC model to a DRN format file.
    #    
    def export_to_drn(self, file):
        stormpy.export_to_drn(self.model, file)



    

