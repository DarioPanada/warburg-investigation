import os
import pandas as pd

from analyzers.SingleReportModelAnalyzers import get_post_execution_analysis
from aws.ExperimentReader import read_experiment_from_queue
from model.models.model_warburg import *

queue_url = "https://sqs.us-east-2.amazonaws.com/746221766782/warburg.fifo"
experiments_file = read_experiment_from_queue(
    queue_url,
    "./experiments",
    num_experiments=1
)

while experiments_file is not None:
    print(experiments_file)

    experiments_file = read_experiment_from_queue(
        queue_url,
        "./experiments",
        num_experiments=1
    )

    experiments = pd.read_csv(experiments_file).to_dict(orient="records")

    num_epochs = 2
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
        print("Uploading to bucket")
        upload_command = "aws s3 sync {0} s3://panaxea-warburg-results/{1}" \
            .format(
            experiment_dir,
            experiment["name"]
        )
        print(upload_command)
        os.system(upload_command)

        rm_command = "rm -r {0}".format(experiment_dir)
        print(rm_command)
        os.system(rm_command)
        print("All done!")
