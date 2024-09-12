import argparse


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

    parser.add_argument(
        "--name",
        type=str,
        default="test",
        help="name for results files"
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=1800,
        help="timeout, in seconds"
    )

    parser.add_argument(
        "--maxiter",
        type=int,
        default=200,
        help="maximum nr of iterations"
    )

    parser.add_argument(
        "--test",
        type=bool,
        default=False,
        help="run tests"
    )

    return parser.parse_args()
