


class Interval:
    def __init__(self, lowerbound, upperbound, state, action, successor, name, id):
        self.lowerbound = lowerbound
        self.upperbound = upperbound
        self.state = state
        self.action = action
        self.successor = successor
        self.name = name
        self.id = id

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



def parse_interval_file(interval_path):
    # adhere to strict naming convention to make parsing the intervals easy
    # could be more efficient, but should be fine for now
    intervals = dict()

    # build intervals from file
    file = open(interval_path,'r')
    for line in file:
        split = line.split()
        interval_name = split[0]
        interval_id = interval_name.split('_')
        state = interval_id[1]
        action = interval_id[2]
        successor = interval_id[3]

        interval = split[1]
        interval = interval.replace('[','')
        interval = interval.replace(']','')
        interval = interval.replace(' ','')
        values = interval.split(',')
        lower = float(values[0])
        upper = float(values[1])

        i = Interval(lower, upper, state, action, successor, interval_name)
        intervals[interval_name] = i
    file.close()

    return intervals


def parse_model_interval(model, pomdp_parameter, interval_path):
    # adhere to strict naming convention to make parsing the intervals easy
    # could be more efficient, but should be fine for now
    intervals = dict()
    items = []
    file = open(interval_path, 'r')
    numline = 0
    for line in file:
        split = line.split()
        interval_name = split[0]

        interval = split[1]
        interval = interval.replace('[', '')
        interval = interval.replace(']', '')
        interval = interval.replace(' ', '')
        values = interval.split(',')
        lower = float(values[0])
        upper = float(values[1])
        i = Item(lower, upper, interval_name)
        items.append(i)


    for state in model.states:
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
                                        i = Interval(item.lowerbound, item.upperbound, str(state.id), 0, str(succ), str(t.monomial[0][0].name), t.monomial[0][0].id)
                                        intervals[str(t.monomial[0][0].name)] = i




    return intervals, items


def parse_model_interval_params(model, intervals_val):
    # adhere to strict naming convention to make parsing the intervals easy
    # could be more efficient, but should be fine for now
    intervals = []

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
                                    if item.name == t.monomial[0][0].name[0]:
                                        i = Interval(item.lowerbound, item.upperbound, int(state.id), 0, int(succ), "I")
                                        intervals.append(i)
                                        break
    return intervals
