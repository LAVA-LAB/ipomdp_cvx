import csv




class ModelInfo:
    def __init__(self):
        self.model_name = ""
        self.nr_model_states = 0
        self.nr_model_actions = 0
        self.nr_model_transitions = 0
        self.mem_size = 0
        self.nr_product_states = 0
        self.nr_product_actions = 0
        self.nr_product_transitions = 0
        self.nr_encoding_variables = 0
        self.nr_encoding_constraints = 0


class Result:

    def __init__(self):
        self.nr_encoding_variables = 0
        self.nr_encoding_constraints = 0
        self.final_result = 0
        self.total_time = 0
        self.steps = []
        self.computation_time = 0
        self.exit_code = ""


    def addStep(self, result, time):
        self.steps.append([result,time])


    def finalize(self):
        self.final_result = self.steps[-1].result
        self.total_time = self.steps[-1].time

    def write_to_files(self, name=""):
        config_name = name+"_config.txt"
        data_name = name+"_results.csv"

        with open(config_name, "w") as f:
            f.write("model name: {0}\n".format(str(self.model_name)))
            f.write("nr model states: {0}\n".format(str(self.nr_model_states)))
            f.write("nr model actions: {0}\n".format(str(self.nr_model_actions)))
            f.write("nr model transitions: {0}\n".format(str(self.nr_model_transitions)))
            f.write("mem size: {0}\n".format(str(self.mem_size)))
            f.write("nr product states: {0}\n".format(str(self.nr_product_states)))
            f.write("nr product actions: {0}\n".format(str(self.nr_product_actions)))
            f.write("nr product transitions: {0}\n".format(str(self.nr_product_transitions)))
            f.write("nr encoding variables: {0}\n".format(str(self.nr_encoding_variables)))
            f.write("nr encoding constraints: {0}\n".format(str(self.nr_encoding_constraints)))
            f.write("final result: {0}\n".format(str(self.final_result)))
            f.write("total time: {0}\n".format(str(self.total_time)))
            f.write("exit code: {0}".format(str(self.exit_code)))


        fields = ["result","time"]
        with open(data_name, 'w') as f:
            write = csv.writer(f)

            write.writerow(fields)
            write.writerows(self.steps)
