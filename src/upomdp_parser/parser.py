import stormpy
import inspect
import interval_parser



def export_to_drn(pomdp, export_file):
    stormpy.export_to_drn(pomdp, export_file)



def parse_prism(prism_file, property_string):
    prism_program = stormpy.parse_prism_program(prism_file)

    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    properties = stormpy.parse_properties_for_prism_program(property_string, prism_program)
    # construct the pPOMDP
    print(inspect.getfullargspec(stormpy.build_parametric_model))
    pomdp = stormpy.build_parametric_model(prism_program, properties)

    pomdp_parameters = pomdp.collect_probability_parameters()

    return pomdp, pomdp_parameters



def parse_drn(drn_file, property_string):
    drn = stormpy.build_parametric_model_from_drn(drn_file)

    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    properties = stormpy.parse_properties_for_prism_program(property_string, drn)
    # construct the pPOMDP
    print(inspect.getfullargspec(stormpy.build_parametric_model))
    pomdp = stormpy.build_parametric_model(drn, properties)

    # get all parameters in the model
    pomdp_parameters = pomdp.collect_probability_parameters()


    return pomdp, pomdp_parameters


def main():
    # basic idea, load the model with parameters on the transitions where intervals should occur,
    # load a second file mapping these parameters to intervals

    upomdp, params = parse_prism("aircraft_small.prism", "Pmax=?[F \"goal\"]")
    intervals, items = interval_parser.parse_model_interval(upomdp, params, "aircraft_small.intervals")

    # loop over the model
    for state in upomdp.states:
        for action in state.actions:
            for transition in action.transitions:
                transition_value = transition.value()
                if transition_value.is_constant():
                    # it's a number!
                    # do whatever you need to do
                    continue
                else:
                    # we assume the denominator is constant, otherwise this breaks
                    poly = transition_value.numerator.polynomial()
                    for term in poly:
                        if term.is_constant():
                            # if the term is a constant, e.g. 0.5 or 1, skip
                            continue
                        else:
                            # we're looking at a variable
                            if not term.monomial[0][0].name == None:
                                name = term.monomial[0][0].name
                                if not str(name) in intervals.keys():
                                    raise RuntimeError("Parameter that was not defined in the interval file")
                                else:
                                    # all looks good, we now have a transition (s,a,s') with the interval given by intervals[term.monomial[0][0].name]
                                    upper_bound = intervals[term.monomial[0][0].name].get_upperbound()
                                    lower_bound = intervals[term.monomial[0][0].name].get_lowerbound()
                                    interval_str = "[{},{}]".format(lower_bound, upper_bound)
                                    print(
                                        f"transition ({state}, {action}, {transition.column}) with parameter {name} has interval {interval_str}")
                            else:
                                continue


if __name__ == "__main__":
    main()