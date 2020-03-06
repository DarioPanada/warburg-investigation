from model.models.model_warburg import *
import pandas as pd

experiments_file = "experiments/experiments_warburg.csv"
experiments = pd.read_csv(experiments_file).to_dict(orient="records")

print("There are {0} experiments".format(len(experiments)))

for experiment in experiments[:5]:
    properties = generate_properties(experiment)
    print(properties)
