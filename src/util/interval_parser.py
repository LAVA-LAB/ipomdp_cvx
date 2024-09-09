import cdd
import numpy as np
import itertools
import stormpy
import stormpy.core
import stormpy.logic
import stormpy.pars
import re
import stormpy.examples
import stormpy.examples.files
import time
import math

from gurobipy import *
# maybe do some type conversion?
class Interval:
    def __init__(self, lowerbound, upperbound, state, action, successor, name):
        self.lowerbound = lowerbound
        self.upperbound = upperbound
        self.state = state
        self.action = action
        self.successor = successor
        self.name = name

    def get_upperbound(self):
        return self.upperbound

    def get_lowerbound(self):
        return self.lowerbound

    def get_state(self):
        return self.state

    def get_action(self):
        return self.action

    def get_successor(self):
        return self.successor

    def __cmp__(self, other):
        return self.state == other.state and self.action == other.action and self.successor == other.action and self.name == other.name

    def __str__(self):
        return "Interval at ({0},{1},{2}) : [{3}, {4}]".format(str(self.state), str(self.action), str(self.successor), str(self.lowerbound), str(self.upperbound))


class Item:
    def __init__(self, lowerbound, upperbound, name):
        self.lowerbound = lowerbound
        self.upperbound = upperbound
        self.name = name

    def get_upperbound(self):
        return self.upperbound

    def get_lowerbound(self):
        return self.lowerbound



 #   def __cmp__(self, other):
 #       return self.state == other.state and self.action == other.action and self.successor == other.action and self.name == other.name

    def __str__(self):
        return "Interval with name {0} with lower and upper bound [{1}, {2}]".format(str(self.name), str(self.lowerbound), str(self.upperbound))


class VertexValue:

    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __str__(self):
        return str(self.value)

class Vertex:
    def __init__(self, vertex_values):
        self.vertex_values = vertex_values

    def get_vertex_value(self, name):
        for v in self.vertex_values:
            if v.name == name:
                return v
        return None

    def __str__(self):
        s = "("
        for p in self.vertex_values:
            s += " {0} ".format(str(p))
        s += ")"
        return s


class Polyhedron:

    def __init__(self, state, action):
        self.state = state
        self.action = action
        self.id = (state, action)
        self.intervals = []
        self.dimension = 0
        self.vertices = []

    def get_state(self):
        return self.state

    def get_action(self):
        return self.action

    def get_id(self):
        return self.id

    def add_interval(self,interval):
        self.intervals.append(interval)
        self.dimension += 1

    def get_dimension(self):
        return self.dimension

    def get_intervals(self):
        return self.intervals

    def contains(self, interval):
        return interval in self.intervals

    def __cmp__(self, other):
        return self.state == other.state and self.action == other.action

    def __str__(self):
        s = "{0} dimensional polyhedron at {1} composed of:\n".format(str(self.dimension), str(self.id))
        for i in self.intervals:
            s += "    {0}\n".format(str(i))
        return s

    def compute_vertices(self):
        # some matrix blocks to help build the inequality matrix
        n = self.dimension
        In = np.identity(n)
        Hn = np.ones((1, n))  # single row matrix of n 1s

        # build matrix A:
        A = np.vstack((-1 * In, In, Hn, -1 * Hn))

        # parse the intervals to build vector b
        ub = []
        lb = []
        for i in self.intervals:
            ub.append([i.upperbound])
            lb.append([i.lowerbound])

        ub = np.array(ub)
        lb = -1 * np.array(lb)
        b = np.vstack((lb, ub, 1, -1))

        # cddlib wants the input Ax <= b as a single matrix [b -A]
        M = np.concatenate((b, -1 * A), axis=1)

        # build cdd matrix, set representation to inequalities
        mat = cdd.Matrix(M, number_type='float')
        mat.rep_type = cdd.RepType.INEQUALITY

        # build polyhedron
        poly = cdd.Polyhedron(mat)

        # get vertices
        gen = poly.get_generators()

        vertices = []
        for i in range(gen.row_size):
            points = []
            #vertex = []
            if gen[i][0] == 0:
                break
            for j in range(gen.col_size - 1):
                value = gen[i][j + 1]
                points.append(VertexValue(self.intervals[j].name, value))
                #vertex.append(value)
            vertex = Vertex(points)
            vertices.append(vertex)

        self.vertices = vertices

        return vertices

    def print_vertices(self):
        s = "["
        for v in self.vertices:
            s += " {0} ".format(str(v))
        s += "]"
        print(s)



def calculate_vertex_combinations(polyhedrons):
    if len(polyhedrons) == 0:
        return []
    elif len(polyhedrons) == 1:
        return [polyhedrons[0].vertices]
    else:
        to_combine = [p.vertices for p in polyhedrons]
        #print(to_combine)
        combinations = list(itertools.product(*to_combine))
        return combinations




def parse_model_interval(model,pomdp_parameter,interval_path):
    # adhere to strict naming convention to make parsing the intervals easy
    # could be more efficient, but should be fine for now
    intervals = []
    polyhedrons = []
    items =[]
    file = open(interval_path, 'r')
    numline = 0
    for line in file:
        split = line.split()
        interval_name = split[0]
        interval_id = interval_name.split('_')
        # state = interval_id[1]
        # action = interval_id[2]
        # successor = interval_id[3]

        interval = split[1]
        interval = interval.replace('[', '')
        interval = interval.replace(']', '')
        interval = interval.replace(' ', '')
        values = interval.split(',')
        lower = float(values[0])
        upper = float(values[1])
        i = Item(lower,upper,interval_name)
        items.append(i)
    for state in model.states:
        #print(state.id,i)
        # print (prob1.get(state))
        # Cons=values constraints on the right hand side for a pdtmc


        for action in state.actions:
            for transition in action.transitions:

                transition_value = transition.value()


                succ = int(transition.column)
                    # Value of transition
                transition_value = transition.value()
                    # TODO check that the denominator is indeed constant?
                    # Denominator of transition
                den = transition_value.denominator.constant_part()
                denom1 = 1 / float(den)

                # If the transition value is not constant
                if not transition_value.is_constant():
                    num = transition_value.numerator.polynomial()

                        # Iterates over terms in numerators
                    for t in num:
                            # If the transition term is a constant
                        if t.is_constant():
                                # Add just value of transition is the successor state is prob1 to the constraints
                            pass
                            # If the transition term is not a constant
                        else:
                            if t.tdeg > 1:
                                raise RuntimeError("We expect the term to be a single variable")
                            if t.monomial[0][0] in pomdp_parameter:
                                for item in items:
                                    if item.name==t.monomial[0][0].name:

                                        i = Interval(item.lowerbound, item.upperbound, int(state.id), 0, int(succ), str(t.monomial[0][0].name))
                                        intervals.append(i)
                                        #sprint(i)
                                        #break
    return intervals, polyhedrons,items

def parse_model(model,intervals_val):
    # adhere to strict naming convention to make parsing the intervals easy
    # could be more efficient, but should be fine for now
    intervals = []
    polyhedrons = []

    for state in model.states:
        # print(state.id)
        # print (prob1.get(state))
        # Cons=values constraints on the right hand side for a pdtmc


        for action in state.actions:
            for transition in action.transitions:

                transition_value = transition.value()


                succ = int(transition.column)
                    # Value of transition
                transition_value = transition.value()
                    # TODO check that the denominator is indeed constant?
                    # Denominator of transition
                den = transition_value.denominator.constant_part()
                denom1 = 1 / float(den)

                # If the transition value is not constant
                if not transition_value.is_constant():
                    num = transition_value.numerator.polynomial()

                        # Iterates over terms in numerators
                    for t in num:
                            # If the transition term is a constant
                        if t.is_constant():
                                # Add just value of transition is the successor state is prob1 to the constraints
                            pass
                            # If the transition term is not a constant
                        else:
                            if t.tdeg > 1:
                                raise RuntimeError("We expect the term to be a single variable")
                            if t.monomial[0][0].name[0] == 'I':
                                for item in intervals_val:
                                    if item.name==t.monomial[0][0].name[0]:
                                        i = Interval(item.lowerbound, item.upperbound, int(state.id), 0, int(succ), "I")
                                        intervals.append(i)
                                        break
    return intervals, polyhedrons



def parse_input(interval_path):
    # adhere to strict naming convention to make parsing the intervals easy
    # could be more efficient, but should be fine for now
    intervals = []
    polyhedrons = []

    # build intervals from file
    file = open(interval_path,'r')
    numline=0
    for line in file:
        split = line.split()
        interval_name = split[0]
        interval_id = interval_name.split('_')
        # state = interval_id[1]
        # action = interval_id[2]
        # successor = interval_id[3]

        interval = split[1]
        interval = interval.replace('[','')
        interval = interval.replace(']','')
        interval = interval.replace(' ','')
        values = interval.split(',')
        lower = float(values[0])
        upper = float(values[1])
        state=41
        action=0
        successor=6
        i = Interval(lower, upper, state, action, successor, interval_name)
        intervals.append(i)
        numline+=1

    i = Interval(lower, upper, 41, 0, 12, "I")
    intervals.append(i)
    i = Interval(lower, upper, 42, 0, 10, "I")
    intervals.append(i)
    i = Interval(lower, upper, 42, 0, 14, "I")
    intervals.append(i)
    file.close()

    # build polyhedrons
    numinterval=0
    for i in intervals:
        numinterval+=1
        if numinterval % 1000==0:
            pass
            #print(numinterval)
        exists = False
        for p in polyhedrons:
            if p.state == i.state and p.action == i.action:
                exists = True
                break
        if not exists:
            polyhedrons.append(Polyhedron(i.state, i.action))

    # put intervals in appropriate polyhedron
    numpolyhedron=0
    for p in polyhedrons:
        numpolyhedron+=1
        if numpolyhedron % 1000==0:
            pass
            #print(numpolyhedron)
        for i in intervals:
            if p.state == i.state and p.action == i.action:
                p.add_interval(i)

    return intervals, polyhedrons



"""
ints, polys = parse_input("/home/marnix/PycharmProjects/Storm/input/pMC1.intervals")
for p in polys:
    p.compute_vertices()
    p.print_vertices()
    for v in p.vertices:
        print(v)
combs = calculate_vertex_combinations(polys)
for inst in combs:
    for v in inst:
        print(type(v))
        for val in v.vertex_values:
            print(type(val))
            print(val.name)
"""
