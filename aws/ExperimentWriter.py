import hashlib
from random import randint

import boto3

sqs = boto3.client('sqs')


def write_experiment_to_queue(queue_url, experiment, exp_name, exp_group,
                              deduplication_id, header,
                              message_group="1"):
    """
    Writes an experiment (provided as a string) to an amazon queue as a
    single message.

    Parameters
    ----------
    queue_url : string
        The url of a queue
    experiment : string
        The csv string representing the experiment
    exp_name : string
        The name of the experiment. Will be provided as an attribute to the
        message sent to the queue. This should be used to indicate the
        "name" property of each csv row.
    exp_group : string
        The group to which the experiment belongs. Will be provided as an
        attribute to the message sent to the queue. This should be the name
        of the csv file to which a row belongs, or another identifier useful
        to re-group experiments if multiple csvs are uploaded to the same
        queue.
    deduplication_id : string
        A unique value, a further random value will be appended to it.
    header : string
        The experiment file's header, will be added as an attribute
    message_group : int, optional
        The group to which the message belongs. Optional, defaults to 1.
    """

    m = hashlib.sha256()
    m.update(str(deduplication_id) + str(randint(0, 100)))
    deduplication_id = m.hexdigest()

    attributes = {
        "experiment_name": {
            "DataType": "String",
            "StringValue": exp_name
        },
        "experiment_group": {
            "DataType": "String",
            "StringValue": exp_group
        },
        "experiment_header": {
            "DataType": "String",
            "StringValue": header
        }
    }

    response = sqs.send_message(
        QueueUrl=queue_url,
        MessageAttributes=attributes,
        MessageBody=experiment,
        MessageGroupId=message_group,
        MessageDeduplicationId=deduplication_id
    )

    print("I have written message {0} with response {1}".format(
        exp_name, response["ResponseMetadata"]["HTTPStatusCode"]))

    if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
        print("Something went wrong! Full response:")
        print(response)


if __name__ == "__main__":

    experiment_file = "experiments_sample_reduced.csv"
    queue_url = "https://sqs.us-east-2.amazonaws.com/746221766782/test-queue" \
                ".fifo"
    group = experiment_file.split(".")[0]
    with open(experiment_file, 'r') as f:
        for n, line in enumerate(f):
            line = line.replace("\n", "")
            if n == 0:
                header = line
            else:
                experiment_name = line.split(",")[0]

                write_experiment_to_queue(queue_url, line, experiment_name,
                                          group, n, header,
                                          message_group=group)
