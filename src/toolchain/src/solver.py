import stormpy
from pomdpdrngurobi_sjabove_robust_interval import QcqpOptions
from pomdpdrngurobi_sjabove_robust_interval import QcqpSolver_affine_simple_fun


def solve_nominal_pomdp(pmc_path, formula_str, threshold, direction="below", maxiter=100, silent=False):
    print("Building model from {}".format(pmc_path))
    pmc = stormpy.build_parametric_model_from_drn(pmc_path)
    fsc_parameters = pmc.collect_probability_parameters()
    print("number of pomdp parameters={}".format(len(fsc_parameters)))
    print("number of pmc states={}".format(pmc.nr_states))

    properties = stormpy.parse_properties(formula_str)

    # compute prob01max states for any policy
    prob0E, prob1A = stormpy.prob01max_states(pmc, properties[0].raw_formula.subformula)
    # solver parameters
    options = QcqpOptions(mu=1e4, maxiter=maxiter, graph_epsilon=1e-2, silent=silent)

    # run uncertain POMDP solver without uncertainty
    solver = QcqpSolver_affine_simple_fun()
    pomdp_parameters = set()
    intervals = {}
    items = {}
    result = solver.run(pmc, fsc_parameters, pomdp_parameters, properties, prob0E, prob1A, threshold, direction,
                        options, intervals, items, True)

    print("SOLVED")
    print("number of iterations={}".format(solver.iterations))
    print("solver time={}".format(solver.solver_timer))
    print("value at initial state={}".format(result.value_at_initial))

    return result


if __name__ == '__main__':
    solve_nominal_pomdp("export_1_mem_collision_partial_obs_2d_upd_hobs_20_small.drn", "P=? [F \"goal\"]", threshold=0.99, silent=False)
