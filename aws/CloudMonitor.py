import boto3
import json
import os
import time

sqs = boto3.client('sqs')


def download_from_queue(queue_url, out_file):
    """
    Downloads up to 5 messages from a queue, outputs their value to a csv in
    format:
    experiment_name, epoch, body
    Parameters
    ----------
    queue_url : string
        The url of the queue
    out_file : string
        Path to the output csv file
    """
    response = sqs.receive_message(
        QueueUrl=queue_url,
        AttributeNames=[
            'experiment_name'
            'epoch',
        ],
        MaxNumberOfMessages=5,
        MessageAttributeNames=[
            'All'
        ],
        VisibilityTimeout=100,
        WaitTimeSeconds=0
    )

    messages = response.get('Messages')

    if messages is None:
        print("No messages found...")
        return

    messages_as_list = []

    for message in messages:
        body = message["Body"]
        attributes = message["MessageAttributes"]
        epoch = attributes["epoch"]["StringValue"]
        experiment_name = attributes["experiment_name"]["StringValue"]

        message_as_list = [experiment_name, epoch, body]
        messages_as_list.append(message_as_list)

        receipt = message['ReceiptHandle']

        sqs.delete_message(
            QueueUrl=queue_url,
            ReceiptHandle=receipt
        )
        print("Experiment: {0}, Epoch: {1}, Body: {2}".format(
            experiment_name,
            epoch,
            body
        ))

    write_header = not os.path.isfile(out_file)

    with open(out_file, "a") as f:
        if write_header:
            header = "experiment_name,epoch,body\n"
            f.write(header)

        for message_as_list in messages_as_list:
            message_string = ",".join(message_as_list) + "\n"
            f.write(message_string)
        f.close()


if __name__ == "__main__":
    with open("../config.json", "r") as f:
        config = json.load(f)
        f.close()

    time_interval = 1

    while True:
        print("Downloading messages...")
        download_from_queue(
            config["aws"]["messages_queue"],
            "cloud_logs/messages_queue.csv",
        )

        print("Downloading exceptions...")
        download_from_queue(
            config["aws"]["exceptions_queue"],
            "cloud_logs/exceptions_queue.csv"
        )

        time.sleep(time_interval)
