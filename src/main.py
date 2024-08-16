from util.argparser import parse_args
import scp.scp_reach_maxmin as scp_prob
import ccp.pomdpdrngurobi_sjabove_robust as ccp_prob
import scp.scp_reward_minmax as scp_rew
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

def run_ccp_prob(model_path, intervals_path, spec, threshold, runname = ""):
    model_info, results = ccp_prob.main(model_path, intervals_path, spec, threshold)
    save_model_info(model_info, "results/info_{0}.txt".format(runname))
    save_results(results, "results/results_{0}.csv".format(runname))

def run_scp_prob_drn(model_path, intervals_path, spec, threshold, memval, runname = ""):
    model_info, results = scp_prob.main_drn(model_path, intervals_path, spec, threshold, memval)
    save_model_info(model_info, "results/info_{0}.txt".format(runname))
    save_results(results, "results/results_{0}.csv".format(runname))


def run_scp_prob_prism(model_path, intervals_path, spec, threshold, memval, runname = ""):
    model_info, results = scp_prob.main_prism(model_path, intervals_path, spec, threshold, memval)
    save_model_info(model_info, "results/info_{0}.txt".format(runname))
    save_results(results, "results/results_{0}.csv".format(runname))


def run_scp_rew_prism(model_path, intervals_path, spec, threshold, memval, runname = ""):
    model_info, results = scp_rew.main_prism(model_path, intervals_path, spec, threshold, memval)
    save_model_info(model_info, "results/info_{0}.txt".format(runname))
    save_results(results, "results/results_{0}.csv".format(runname))


def run_scp_rew_drn(model_path, intervals_path, spec, threshold, memval, runname = ""):
    model_info, results = scp_rew.main_drn(model_path, intervals_path, spec, threshold, memval)
    save_model_info(model_info, "results/info_{0}.txt".format(runname))
    save_results(results, "results/results_{0}.csv".format(runname))


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
    runname = args['runname']



    if method == "ccp":
        if drn == "":
            print("CCP method only supports drn files")
            exit(0)
        model_path = drn
        if tp == "prob":
            run_ccp_prob(model_path, intervals_path, spec, threshold)
        elif tp == "rew":
            print("reward solver for CCP method not supported")
            exit(0)
        else:
            print("Unknown or incorrect property type given, ensure --type is either \"prob\" or \"rew\"")
            exit(0)
    else:
        if tp == "prob":
            if drn == "":
                run_scp_prob_prism(prism, intervals_path, spec, threshold, memval, runname)
            else:
                run_scp_prob_drn(drn, intervals_path, spec, threshold, memval, runname)
        elif tp == "rew":
            if drn == "":
                run_scp_rew_prism(prism, intervals_path, spec, threshold, memval, runname)
            else:
                run_scp_rew_drn(drn, intervals_path, spec, threshold, memval, runname)
        else:
            print("Unknown or incorrect property type given, ensure --type is either \"prob\" or \"rew\"")
            exit(0)
