import boto3
import datetime
import random
import time

from Common import get_instance_and_spot_request_id


def write_message_to_queue(queue_url, experiment_name, message_text,
                           model=None):
    """
    Writes an experiment to an aws queue.

    Parameters
    ----------
    queue_url : string
        AWS queue url
    experiment_name : string
        Name of the experiment that raised the exception
    message_text : string
        Message to write
    model : Model, optional
        Current model instance when the exception was raised. Optional,
        default to none
    """
    sqs = boto3.client('sqs')

    instance_id, request_id = get_instance_and_spot_request_id()

    if model is not None:
        epoch = model.current_epoch
    else:
        epoch = None

    attributes = {
        "experiment_name": {
            "DataType": "String",
            "StringValue": experiment_name
        },
        "epoch": {
            "DataType": "String",
            "StringValue": str(epoch)
        },
        "timestamp": {
            "DataType": "String",
            "StringValue": str(datetime.datetime.now().isoformat())
        },
        "instance_id": {
            "DataType": "String",
            "StringValue": str(instance_id)
        },
        "request_id": {
            "DataType": "String",
            "StringValue": str(request_id)
        }
    }

    # The message group id is the general name of the experiment, minus the
    # suffix indicating the simulation number
    message_group_id = "_".join(experiment_name.split("_")[:-1])

    response = sqs.send_message(
        QueueUrl=queue_url,
        MessageAttributes=attributes,
        MessageBody=message_text,
        MessageGroupId=message_group_id,
        MessageDeduplicationId=str(time.time()) + str(random.random())

    )

    print("I have written message {0} with response {1}".format(
        message_text, response["ResponseMetadata"]["HTTPStatusCode"]))

    if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
        print("Something went wrong! Full response:")
        print(response)
