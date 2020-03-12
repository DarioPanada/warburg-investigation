import boto3
import time
import random
sqs = boto3.client('sqs')


def write_message_to_queue(queue_url, experiment_name, exception_text,
                           model):
    """
    Writes an experiment to an aws queue.

    Parameters
    ----------
    queue_url : string
        AWS queue url
    experiment_name : string
        Name of the experiment that raised the exception
    exception_text : string
        String representation of the exception
    model : Model
        Current model instance when the exception was raised.
    """

    attributes = {
        "experiment_name": {
            "DataType": "String",
            "StringValue": experiment_name
        },
        "epoch": {
            "DataType": "String",
            "StringValue": str(model.current_epoch)
        }
    }

    # The message group id is the general name of the experiment, minus the
    # suffix indicating the simulation number
    message_group_id = "_".join(experiment_name.split("_")[:-1])

    response = sqs.send_message(
        QueueUrl=queue_url,
        MessageAttributes=attributes,
        MessageBody=exception_text,
        MessageGroupId=message_group_id,
        MessageDeduplicationId=str(time.time())+str(random.random())

    )

    print("I have written message {0} with response {1}".format(
        exception_text, response["ResponseMetadata"]["HTTPStatusCode"]))

    if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
        print("Something went wrong! Full response:")
        print(response)
