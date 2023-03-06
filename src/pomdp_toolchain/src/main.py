from pomdp_solver import *
import stormpy.examples
import stormpy.examples.files
import stormpy.pars
from pycarl.cln.cln import Rational
import time

def prism_maze():
    path = "../models/maze.prism"
    prop ="P=?[F \"goal\" ]"
    policy = solve_pomdp(path, prop)
    for value in policy.parameter_values:
        print("{}: {}".format(value, policy.parameter_values[value]))

def drn_maze():
    print("***** Starting main_drn *****")
    path = "../models/maze.drn"
    prop = "P=?[F \"goal\"]"
    policy = solve_pomdp(path, prop)

def prism_maze2():
    path = "../models/maze2.prism"
    prop = "P=?[F \"target\"]"
    policy = solve_pomdp(path, prop)

def drn_maze2():
    print("***** Starting main_drn *****")
    path = "../models/maze2.drn"
    prop = "P=?[F \"target\"]"
    policy = solve_pomdp(path, prop)

def prism_3x3grid():
    path = "../models/3x3grid.prism"
    prop = "Rmin=? [ F \"goal\" ]"
    policy = solve_pomdp(path, prop)

def drn_3x3grid():
    path = "../models/3x3grid.drn"
    prop = "P=? [ F \"goal\" ]"
    policy = solve_pomdp(path, prop)

def prism_4x4grid():
    path = "../models/4x4grid.prism"
    prop = "Rmin=? [ F \"goal\" ]"
    policy = solve_pomdp(path, prop)

def drn_4x4grid():
    path = "../models/4x4grid.drn"
    prop = "P=? [ F \"goal\" ]"
    policy = solve_pomdp(path, prop)

def prism_repudiation():
    path = "../models/repudiation_pomdp.prism"
    prop = "Pmax=? [ F \"unfair\" ]"
    policy = solve_pomdp(path, prop, undefined_constants="K=4")

def prism_aircraft_small():
    path = "../models/aircraft_small.prism"
    prop = "P=? [F \"goal\"]"
    policy = solve_pomdp(path, prop)

def drn_aircraft_small():
    path = "../models/aircraft_small.drn"
    prop = "P=? [F \"goal\"]"
    policy = solve_pomdp(path, prop)

def prism_aircraft():
    path = "../models/aircraft.prism"
    prop = "P=? [F \"goal\"]"
    policy = solve_pomdp(path, prop)

def drn_aircraft():
    path = "../models/aircraft.prism"
    prop = "P=? [F \"goal\"]"
    policy = solve_pomdp(path, prop)

if __name__ == "__main__":
    #prism_maze()
    #drn_maze()
    #prism_maze2()
    #drn_maze2()
    #prism_3x3grid()
    #drn_3x3grid()
    #prism_4x4grid()
    #drn_4x4grid()
    #prism_repudiation()
    #prism_aircraft_small()
    #drn_aircraft_small()
    #prism_aircraft()
    #drn_aircraft()


    
