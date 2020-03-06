from model.models.model_warburg import *
import pandas as pd

experiments_file = "experiments/experiments_warburg.csv"
experiments = pd.read_csv(experiments_file).to_dict(orient="records")

num_epcohs = 100

print("There are {0} experiments".format(len(experiments)))

for experiment in experiments[:5]:

    print("Running {0}".format(experiment["name"]))
    print("Configuring...")
    properties = generate_properties(experiment)
    print("Generating model...")
    model = generate_model(properties, 5)
    print("Simulating...")
    model.run()
    print("All done!")
