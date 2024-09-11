from util.argparser import parse_args
import scp.scp_reach_maxmin as scp_prob
import ccp.ccp_reach_maxmin as ccp_prob
import scp.scp_reward_minmax as scp_rew
import ccp.ccp_reward_minmax as ccp_rew
import csv



def save_model_info(model_info, filename=""):

    with open(filename, "w") as f:
        for key, val in model_info.items():
            f.write("{0}:    {1}\n".format(str(key), str(val)))


def save_results(results, filename=""):
    fields = ["result", "time"]
    with open(filename, 'w') as f:
        write = csv.writer(f)
        write.writerow(fields)
        write.writerows(results)


def save_eval_results(eval_results, filename=""):

    with open(filename, "w") as f:
        for key, val in eval_results.items():
            f.write("{0}:    {1}\n".format(str(key), str(val)))




def run_ccp_prob_prism(model_path, intervals_path, spec, threshold, memval, runname = "", timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = ccp_prob.main_prism(model_path, intervals_path, spec, threshold, memval, timeout, maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))


def run_ccp_prob_drn(model_path, intervals_path, spec, threshold, memval, runname = "", timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = ccp_prob.main_drn(model_path, intervals_path, spec, threshold, memval, timeout, maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))


def run_ccp_rew_prism(model_path, intervals_path, spec, threshold, memval, runname = "", timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = ccp_rew.main_prism(model_path, intervals_path, spec, threshold, memval, timeout, maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))


def run_ccp_rew_drn(model_path, intervals_path, spec, threshold, memval, runname = "", timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = ccp_rew.main_drn(model_path, intervals_path, spec, threshold, memval, timeout, maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))

def run_scp_prob_drn(model_path, intervals_path, spec, threshold, memval, runname = "",timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = scp_prob.main_drn(model_path, intervals_path, spec, threshold, memval,timeout,maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))


def run_scp_prob_prism(model_path, intervals_path, spec, threshold, memval, runname = "",timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = scp_prob.main_prism(model_path, intervals_path, spec, threshold, memval,timeout,maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))


def run_scp_rew_prism(model_path, intervals_path, spec, threshold, memval, runname = "",timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = scp_rew.main_prism(model_path, intervals_path, spec, threshold, memval,timeout,maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))


def run_scp_rew_drn(model_path, intervals_path, spec, threshold, memval, runname = "",timeout=1800, maxiter=200, evaluation_set = []):
    model_info, results, eval_results = scp_rew.main_drn(model_path, intervals_path, spec, threshold, memval,timeout, maxiter, evaluation_set)
    save_model_info(model_info, "results/run_info_{0}.txt".format(runname))
    save_results(results, "results/intermediate_results_{0}.csv".format(runname))
    save_eval_results(eval_results, "results/evaluation_results_{0}.txt".format(runname))




def test_all():
    print("testing all methods")
    print("--------------------\n\n\ntesting scp prob\n\n\n")
    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals","Pmax=?[F \"goal\"]",
                     0.99 ,1, "test_scp_prob",
                     timeout=1800, maxiter=10,
                     evaluation_set=["models/aircraft/aircraft_small.intervals", "models/aircraft/aircraft_nominal.intervals"])


    print("--------------------\n\n\ntesting ccp prob\n\n\n")
    run_ccp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 1, "test_ccp_prob",
                     timeout=1800, maxiter=10,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_nominal.intervals"])

    print("--------------------\n\n\ntesting scp rew\n\n\n")
    run_scp_rew_prism("models/intercept/intercept3.prism", "models/intercept/intercept_big.intervals", "Rmin=?[F \"goal\"]",
                      threshold=2, memval=1,runname="test_scp_rew",timeout=1800, maxiter=10,
                      evaluation_set=["models/intercept/intercept_small.intervals","models/intercept/intercept_nominal.intervals"])

    print("--------------------\n\n\ntesting ccp rew\n\n\n")
    run_ccp_rew_prism("models/intercept/intercept3.prism", "models/intercept/intercept_big.intervals", "Rmin=?[F \"goal\"]",
                      threshold=2, memval=1,runname="test_ccp_rew",timeout=1800, maxiter=10,
                      evaluation_set=["models/intercept/intercept_small.intervals","models/intercept/intercept_nominal.intervals"])

    print("--------------------\n\n\nDONE")


if __name__ == '__main__':
    args = vars(parse_args())
    prism = args['prism']
    drn = args['drn']
    intervals_path = args['intervals']
    tp = args['type']
    spec = args['spec']
    threshold = args['threshold']
    method = args['method']
    memval = args['memval']
    runname = args['name']
    timeout = args['timeout']
    maxiter = args['maxiter']
    test = args['test']

    if test:
        test_all()
        exit()


    if method == "ccp":
        if tp == "prob":
            if drn == "":
                run_ccp_prob_prism(prism, intervals_path, spec, threshold, memval, runname, timeout, maxiter)
            else:
                run_ccp_prob_drn(drn, intervals_path, spec, threshold, memval, runname, timeout, maxiter)
        elif tp == "rew":
            if drn == "":
                run_ccp_rew_prism(prism, intervals_path, spec, threshold, memval, runname, timeout, maxiter)
            else:
                run_ccp_rew_drn(drn, intervals_path, spec, threshold, memval, runname, timeout, maxiter)
        else:
            print("Unknown or incorrect property type given, ensure --type is either \"prob\" or \"rew\"")
            exit(0)
    else:
        if tp == "prob":
            if drn == "":
                run_scp_prob_prism(prism, intervals_path, spec, threshold, memval, runname, timeout,maxiter)
            else:
                run_scp_prob_drn(drn, intervals_path, spec, threshold, memval, runname, timeout,maxiter)
        elif tp == "rew":
            if drn == "":
                run_scp_rew_prism(prism, intervals_path, spec, threshold, memval, runname, timeout,maxiter)
            else:
                run_scp_rew_drn(drn, intervals_path, spec, threshold, memval, runname, timeout,maxiter)
        else:
            print("Unknown or incorrect property type given, ensure --type is either \"prob\" or \"rew\"")
            exit(0)
