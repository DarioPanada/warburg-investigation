import time
from collections import defaultdict

import boto3
import pandas as pd

sqs = boto3.client('sqs')


def read_experiment_from_queue(queue_url, experiments_dir, num_experiments=1):
    """
    Reads experiments from queue, writes them to csv files
    Parameters
    ----------
    queue_url : string
        The url of the aws queue
    experiments_dir : string
        The directory where the csv files will be saved
    num_experiments : int
        Number of messages to download
    """
    # Retrieving experiments
    response = sqs.receive_message(
        QueueUrl=queue_url,
        AttributeNames=[
            'experiment_name',
            'experiment_group',
            'experiment_header'
        ],
        MaxNumberOfMessages=num_experiments,
        MessageAttributeNames=[
            'All'
        ],
        VisibilityTimeout=100,
        WaitTimeSeconds=0
    )
    messages = response.get('Messages')

    if messages is None:
        print("No messages on queue")
        return None

    print("I have retrieved {0} messages".format(len(messages)))

    experiment_groups = defaultdict(list)

    # Getting important information from the experiments
    for message in messages:
        experiment = {
            "experiment": message["Body"],
            "name": message["MessageAttributes"]["experiment_name"][
                "StringValue"],
            "experiment_group": message["MessageAttributes"][
                "experiment_group"]["StringValue"],
            "experiment_header": message["MessageAttributes"][
                "experiment_header"]["StringValue"]
        }

        experiment_groups[experiment["experiment_group"]].append(experiment)

    experiment_group_names = experiment_groups.keys()

    print("There are {0} experiment groups in this batch, each will be "
          "written to a separate file".format(len(experiment_group_names)))

    for egn in experiment_group_names:
        experiments = experiment_groups[egn]
        # Writing experiment to file
        experiment_group = experiments[0]["experiment_group"]
        experiment_header = experiments[0]["experiment_header"].split(",")

        experiments_row = [e["experiment"].split(",") for e in experiments]

        experiments_df = pd.DataFrame(experiments_row,
                                      columns=experiment_header)

        output_path = "{0}/{1}_{2}.csv".format(experiments_dir,
                                               experiment_group, time.time())

        print("Writing to {0}".format(output_path))
        experiments_df.to_csv(output_path, index=False)

    print("Deleting from queue")

    for message in messages:
        receipt = message['ReceiptHandle']

        sqs.delete_message(
            QueueUrl=queue_url,
            ReceiptHandle=receipt
        )

    print("...all done")


if __name__ == "__main__":
    queue_url = "https://sqs.us-east-2.amazonaws.com/746221766782/test-queue" \
                ".fifo"
    read_experiment_from_queue(queue_url, ".", num_experiments=1)
