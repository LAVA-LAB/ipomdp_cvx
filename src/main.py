import argparse
import scp.pomdpdrngurobi_sjabove_robust_interval as scp_prob
import ccp.pomdpdrngurobi_sjabove_robust as ccp_prob
import scp.pomdpdrngurobi_sjabove_robust_interval_rew as scp_rew

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--prism",
        type=str,
        default="",
        help="path to prism file containing the (u)POMDP"
    )

    parser.add_argument(
        "--drn",
        type=str,
        default="",
        help="path to drn file containing the (u)POMDP"
    )

    parser.add_argument(
        "--intervals",
        type=str,
        help="path to file containing the intervals"
    )

    parser.add_argument(
        "--type",
        type=str,
        default="prob",
        help="prob or rew property"
    )

    parser.add_argument(
        "--spec",
        type=str,
        help="specification, either of the form R=?(F T), or P=?(F T)"
    )

    parser.add_argument(
        "--threshold",
        type=float,
        help="specification threshold, as float"
    )

    parser.add_argument(
        "--method",
        type=str,
        default="scp",
        help="scp or ccp mode"
    )

    parser.add_argument(
        "--memval",
        type=int,
        default=1,
        help="memory size of the policy (SCP method only)"
    )

    return parser.parse_args()


def run_ccp_prob(model_path, intervals_path, spec, threshold):
    ccp_prob.main(model_path, intervals_path, spec, threshold)


def run_scp_prob_drn(model_path, intervals_path, spec, threshold, memval):
    scp_prob.main_drn(model_path, intervals_path, spec, threshold, memval)


def run_scp_prob_prism(model_path, intervals_path, spec, threshold, memval):
    scp_prob.main_prism(model_path, intervals_path, spec, threshold, memval)


def run_scp_rew_prism(model_path, intervals_path, spec, threshold, memval):
    scp_rew.main_prism(model_path, intervals_path, spec, threshold, memval)


def run_scp_rew_drn(model_path, intervals_path, spec, threshold, memval):
    scp_rew.main_drn(model_path, intervals_path, spec, threshold, memval)


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
                run_scp_prob_prism(prism, intervals_path, spec, threshold, memval)
            else:
                run_scp_prob_drn(drn, intervals_path, spec, threshold, memval)
        elif tp == "rew":
            if drn == "":
                run_scp_rew_prism(prism, intervals_path, spec, threshold, memval)
            else:
                run_scp_rew_drn(drn, intervals_path, spec, threshold, memval)
        else:
            print("Unknown or incorrect property type given, ensure --type is either \"prob\" or \"rew\"")
            exit(0)
