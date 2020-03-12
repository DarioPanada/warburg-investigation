import os
import pandas as pd
from analyzers.SingleReportModelAnalyzers import get_post_execution_analysis
from model.models.model_warburg import *
import json

with open("config.json", "r") as f:
    config = json.load(f)
    f.close()

experiments_file = "{0}/{1}".format(
    config["experiments_dir"],
    config["experiment_file"])

experiments = pd.read_csv(experiments_file).to_dict(orient="records")

num_epochs = config["num_epochs"]
output_dir = config["output_dir"]

print("There are {0} experiments".format(len(experiments)))

for experiment in experiments:

    experiment_dir = "{0}/{1}".format(output_dir, experiment["name"])

    if not os.path.isdir(experiment_dir):
        os.mkdir(experiment_dir)

    print("Running {0}".format(experiment["name"]))
    print("Configuring...")
    properties = generate_properties(experiment)
    properties["outDir"] = experiment_dir
    properties["config"] = config
    print("Generating model...")
    model = generate_model(properties, num_epochs)
    print("Simulating...")
    model.run()
    print("Running analysis...")
    get_post_execution_analysis(experiment_dir)

    print("All done!")
