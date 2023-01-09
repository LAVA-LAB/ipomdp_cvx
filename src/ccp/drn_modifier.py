import stormpy
import stormpy.core
import stormpy.logic
import stormpy.pars
import stormpy.pomdp
import re


def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def prism_to_drn(prism_file):
    prism_program = stormpy.parse_prism_program(prism_file)

    formula_str = "P=? [F \"goal\"]"
    interval_path = "aircraft.intervals"
    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = False
    options = stormpy.DirectEncodingOptions()
    options.allow_placeholders = False

    properties = stormpy.parse_properties_for_prism_program(formula_str, prism_program)
    # construct the pPOMDP
    import inspect
    #print(inspect.getfullargspec(stormpy.build_parametric_model))
    pomdp = stormpy.build_parametric_model(prism_program, properties)

    pomdp_parameters = pomdp.collect_probability_parameters()
    stormpy.export_parametric_to_drn(pomdp, "collision_partial_obs_2d_upd_hobs_20_small_2.drn", options)


def parse_model(input_drn="", input_intervals="", output_drn="", output_intervals=""):
    model_path = "export_1_mem_ppomdp.drn"
    interval_path = "satellite.intervals"
    formula_str = "P=? [F \"goal\"]"

    opts = stormpy.DirectEncodingParserOptions()
    opts.build_choice_labels = True
    pomdp_drn_model = stormpy.build_parametric_model_from_drn(model_path, opts)

    interval_parameters = pomdp_drn_model.collect_probability_parameters()

    # print(interval_parameters)

    #pomdp = stormpy.pomdp.make_canonic(pomdp_drn_model)
    pomdp = pomdp_drn_model
    replace = {}
    sets = set()

    for s in pomdp.states:
        for a in s.actions:
            for t in a.transitions:
                fun = t.value()
                num = fun.numerator
                den = fun.denominator
                if not num.is_constant():
                    poly = num.polynomial()
                    if poly.nr_terms > 1:
                        print(type(poly))
                        vars = poly.gather_variables()
                        if vars not in sets:
                            sets.add(frozenset(vars))

    return sets


def replace_param(sets, param_string):
    return False


def replacer(string):
    repl = string.replace("-", "_")
    repl = repl.replace("(", "_")
    repl = repl.replace(")", "_")
    repl = repl.replace("/", "_")
    repl = repl.replace("+", "_")
    repl = repl.replace("*", "_")
    repl = repl.replace(" ", "")
    return repl


def modify_intervals(new_parameters, input_intervals, output_intervals):
    file = open(input_intervals)
    lines = file.readlines()
    file.close()
    output = ""
    old_parameter_mapping = {}
    old_parameters = set()

    new_parameter_mapping = ""

    for line in lines:
        output += line
        split = line.split(" ")
        param = split[0]
        old_parameters.add(param)
        interval = split[1].replace("[","")
        interval = interval.replace("]","")
        bounds = interval.split(",")
        lb = float(bounds[0])
        ub = float(bounds[1])
        old_parameter_mapping[param] = (lb,ub)

    for new_param in new_parameters:
        bound1 = 1
        bound2 = 1

        for old_param in old_parameters:
            if old_param in new_param:
                bound1 -= old_parameter_mapping[old_param][0]
                bound2 -= old_parameter_mapping[old_param][1]

        lb = min(bound1, bound2)
        ub = max(bound1, bound2)
        new_parameter_mapping += new_param.replace("\n"," ") + "[{0},{1}]\n".format(str(lb),str(ub))

    output += new_parameter_mapping

    outputf = open(output_intervals, "w")
    outputf.write(output)
    outputf.close()


def modify_drn(input_drn, output_drn, input_intervals, output_intervals):
    #sets = parse_model(input_drn)

    output = ""
    parameters_string = ""
    model = False
    new_parameters = []

    with open(input_drn) as inp:
        line = inp.readline()

        while(line):
            if not model:
                if line.startswith("@parameters"):
                    # process current parameters
                    output += "@parameters\n"
                    line = inp.readline()
                    parameters = line
                    parameters = parameters.replace("\n", "")
                    parameters_string = parameters
                    output += "$$PARAMETERS$$\n"
                    line = inp.readline()
                elif line.startswith("@model"):
                    model = True
                    output += line
                    line = inp.readline()
                else:
                    output += line
                    line = inp.readline()
            else:
                if "-" in line and not "{" in line:
                    splitted = line.split(":")
                    repl = replacer(splitted[1])
                    param = "I_"+repl
                    if param not in new_parameters:
                        new_parameters.append(param)
                        parameters_string = param + " " + parameters_string
                    new = splitted[0] + ": " + param
                    output += new
                    line = inp.readline()
                else:
                    output += line
                    line = inp.readline()
        if "-" in line and not "{" in line:
            splitted = line.split(":")
            repl = replacer(splitted[1])
            param = "I_" + repl
            parameters = param + " " + parameters
            new = splitted[0] + " : " + param
            output += new
        else:
            output += line
    parameters_string = parameters_string.replace("\n","")
    output = output.replace("$$PARAMETERS$$",parameters_string)
    outputf = open(output_drn, "w")
    outputf.write(output)
    outputf.close()

    modify_intervals(new_parameters,input_intervals,output_intervals)

    print("DONE")





prism_to_drn("collision_partial_obs_2d_upd_hobs_20_small.prism")
#modify_drn("ppomdp.drn","ppomdp_auto_mod.drn","satellite.intervals","satellite_auto_mod.intervals")
















