from main import run_scp_rew_prism
from main import run_scp_rew_drn
from main import run_scp_prob_drn
from main import run_scp_prob_prism
from main import run_ccp_prob_drn
from main import run_ccp_rew_drn
from main import run_ccp_rew_prism



def ex2_intercept_nominal():
    run_scp_rew_prism("models/intercept/intercept3.prism", "models/intercept/intercept_nominal.intervals",
                      "Rmin=?[F \"goal\"]", threshold=2, memval=1, runname="ex2_scp_intercept_nominal_mem1", timeout=1800,
                      maxiter=100, evaluation_set=["models/intercept/intercept_small.intervals",
                                                   "models/intercept/intercept_big.intervals"])

    run_scp_rew_prism("models/intercept/intercept3.prism", "models/intercept/intercept_nominal.intervals",
                      "Rmin=?[F \"goal\"]", threshold=2, memval=2, runname="ex2_scp_intercept_nominal_mem2", timeout=1800,
                      maxiter=100, evaluation_set=["models/intercept/intercept_small.intervals",
                                                   "models/intercept/intercept_big.intervals"])


def ex2_eps_aircraft_nominal():
    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_nominal_eps.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 1, "ex2_scp_aircraft_mem1_eps", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_big.intervals"])

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_nominal_eps.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 2, "ex2_scp_aircraft_mem2_eps", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_big.intervals"])

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_nominal_eps.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 3, "ex2_scp_aircraft_mem3_eps", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_big.intervals"])


def ex2_aircraft_nominal():
    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_nominal.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 1, "ex2_scp_aircraft_mem1", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_big.intervals"])

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_nominal.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 2, "ex2_scp_aircraft_mem2", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_big.intervals"])

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_nominal.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 3, "ex2_scp_aircraft_mem3", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_big.intervals"])




def ex2_eps_satellite_nominal():

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_nominal_eps.intervals", "Rmin=? [F \"goal\"]", 1, 1, "ex2_scp_satellite_mem1",1800,100,["models/satellite/satellite_small.intervals","models/satellite/satellite_big.intervals"])

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_nominal_eps.intervals", "Rmin=? [F \"goal\"]",
                    1, 2, "ex2_scp_satellite_mem2", 1800, 100,
                    ["models/satellite/satellite_small.intervals", "models/satellite/satellite_big.intervals"])

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_nominal_eps.intervals", "Rmin=? [F \"goal\"]", 1, 3, "ex2_scp_satellite_mem3",1800,100,["models/satellite/satellite_small.intervals","models/satellite/satellite_big.intervals"])

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_nominal_eps.intervals", "Rmin=? [F \"goal\"]",
                    1, 4, "ex2_scp_satellite_mem4", 1800, 100,
                    ["models/satellite/satellite_small.intervals", "models/satellite/satellite_big.intervals"])



def ex1_eps_satellite():
    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_big.intervals", "Rmin=? [F \"goal\"]", 1, 1, "ex1_scp_satellite_mem1_eps",1800,100,["models/satellite/satellite_small.intervals","models/satellite/satellite_nominal_eps.intervals"])

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_big.intervals", "Rmin=? [F \"goal\"]",
                    1, 2, "ex1_scp_satellite_mem2_eps", 1800, 100,
                    ["models/satellite/satellite_small.intervals", "models/satellite/satellite_nominal_eps.intervals"])

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_big.intervals", "Rmin=? [F \"goal\"]", 1, 3, "ex1_scp_satellite_mem3_eps",1800,100,["models/satellite/satellite_small.intervals","models/satellite/satellite_nominal_eps.intervals"])

    run_scp_rew_drn("models/satellite/satellite.drn", "models/satellite/satellite_big.intervals", "Rmin=? [F \"goal\"]",
                    1, 4, "ex1_scp_satellite_mem4_eps", 1800, 100,
                    ["models/satellite/satellite_small.intervals", "models/satellite/satellite_nominal_eps.intervals"])



def ex1_intercept():
    run_scp_rew_prism("models/intercept/intercept3.prism", "models/intercept/intercept_big.intervals", "Rmin=?[F \"goal\"]", threshold=2, memval=1,runname="ex1_scp_intercept_mem0",timeout=1800, maxiter=100, evaluation_set=["models/intercept/intercept_small.intervals","models/intercept/intercept_nominal.intervals"])

    print("\n\n----------------\n\n")

    run_ccp_rew_prism("models/intercept/intercept3.prism", "models/intercept/intercept_big.intervals",
                      "Rmin=?[F \"goal\"]", threshold=2, memval=1, runname="ex1_ccp_intercept_mem0", timeout=1800,
                      maxiter=100, evaluation_set=["models/intercept/intercept_small.intervals",

                                                  "models/intercept/intercept_nominal.intervals"])
def ex1_aircraft():

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals","Pmax=?[F \"goal\"]",0.99 ,1, "ex1_scp_aircraft_mem1", timeout=1800, maxiter=100, evaluation_set=["models/aircraft/aircraft_small.intervals", "models/aircraft/aircraft_nominal.intervals"])

    print("\n\n----------------\n\n")

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 2, "ex1_scp_aircraft_mem2", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_nominal.intervals"])

    print("\n\n----------------\n\n")

    run_scp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 3, "ex1_scp_aircraft_mem3", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_nominal.intervals"])

    print("\n\n----------------\n\n")

    run_ccp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 1, "ex1_ccp_aircraft_mem1", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_nominal.intervals"])

    print("\n\n----------------\n\n")

    run_ccp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 2, "ex1_ccp_aircraft_mem2", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_nominal.intervals"])

    print("\n\n----------------\n\n")

    run_ccp_prob_drn("models/aircraft/aircraft.drn", "models/aircraft/aircraft_big.intervals", "Pmax=?[F \"goal\"]",
                     0.99, 3, "ex1_ccp_aircraft_mem3", timeout=1800, maxiter=100,
                     evaluation_set=["models/aircraft/aircraft_small.intervals",
                                     "models/aircraft/aircraft_nominal.intervals"])


def all():
    ex1_intercept()
    #ex1_aircraft()
    ex1_eps_satellite()
    ex2_eps_satellite_nominal()
    #ex2_eps_aircraft_nominal()
    #ex2_aircraft_nominal()
    #ex2_intercept_nominal()

if __name__ == "__main__":
    all()

