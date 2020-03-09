from model.models.model_warburg import *
from analyzers.SingleReportModelAnalyzers import get_post_execution_analysis
import pandas as pd
import os

experiments_file = "experiments/experiments_warburg.csv"
experiments = pd.read_csv(experiments_file).to_dict(orient="records")

num_epochs = 1
output_dir = "reports"

print("There are {0} experiments".format(len(experiments)))

for experiment in experiments[:5]:

    experiment_dir = "{0}/{1}".format(output_dir, experiment["name"])

    if not os.path.isdir(experiment_dir):
        os.mkdir(experiment_dir)

    print("Running {0}".format(experiment["name"]))
    print("Configuring...")
    properties = generate_properties(experiment)
    properties["outDir"] = experiment_dir
    print("Generating model...")
    model = generate_model(properties, num_epochs)
    print("Simulating...")
    model.run()
    print("Running analysis...")
    get_post_execution_analysis(experiment_dir)
    print("All done!")

