import json
import os
import pandas as pd
import time

from analyzers.SingleReportModelAnalyzers import get_post_execution_analysis
from aws.Common import terminate_instance_and_spot_request, \
    get_instance_and_spot_request_id
from aws.ExperimentReader import read_experiment_from_queue
from aws.MessageWriter import write_message_to_queue
from model.models.model_warburg import *

with open("config.json", "r") as f:
    config = json.load(f)
    f.close()

queue_url = config["aws"]["experiments_queue"]
experiments_dir = config["experiments_dir"]
output_dir = config["output_dir"]
num_epochs = config["num_epochs"]

run_analysis = True
retry = 0

while True:

    experiments_file = read_experiment_from_queue(
        queue_url,
        experiments_dir,
        num_experiments=1
    )

    max_retry = 3
    retry_interval = 5

    # If there are no further experiments on the queue
    if experiments_file is None:
        num_active_experiments = len([f for f in os.listdir(output_dir)
                                      if "gitkeep" not in f])
        # And no other process is working on an experiment
        if num_active_experiments == 0:
            # Terminate the instance and the associated spot request
            if retry == max_retry:
                instance_id, request_id = get_instance_and_spot_request_id()
                terminate_instance_and_spot_request(instance_id, request_id)
            else:
                retry += 1
                time.sleep(retry_interval)
    else:
        retry = 0
        experiments_from_queue = pd.read_csv(experiments_file).to_dict(
            orient="records")

        print("There are {0} experiments".format(len(experiments_from_queue)))

        for experiment in experiments_from_queue:
            try:
                experiment_dir = "{0}/{1}".format(output_dir,
                                                  experiment["name"])

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

                write_message_to_queue(
                    config["aws"]["messages_queue"],
                    model.properties["name"],
                    "COMPLETED",
                    model
                )

                if run_analysis:
                    print("Running analysis...")
                    get_post_execution_analysis(experiment_dir)
                print("Uploading to bucket")
                upload_command = "aws s3 sync {0} " \
                                 "s3://panaxea-warburg-results/{1}" \
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
            except Exception as e:
                write_message_to_queue(
                    config["aws"]["exceptions_queue"],
                    experiment["name"],
                    str(e),
                    model
                )
                write_message_to_queue(
                    config["aws"]["messages_queue"],
                    model.properties["name"],
                    "EXCEPTION",
                    model
                )
                rm_command_report = "rm -r {0}".format(experiment_dir)
                print(rm_command_report)
                os.system(rm_command_report)
                print("All done...!")

        rm_command_exp_file = "rm {0}".format(experiments_file)
        print(rm_command_exp_file)
        os.system(rm_command_exp_file)
