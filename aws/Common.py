import boto3
import json
import os

from MessageWriter import write_message_to_queue


def get_instance_and_spot_request_id():
    """
    If running on an aws instance, returns the instance id and spot request id.
    Returns
    -------
    instance_id : string
        The instance id
    request_id : string
        The spot request id

    Or None, None it the application is not running on an aws instance
    """
    user = os.environ.get("user")

    if user == "ec2-user":
        ec2 = boto3.client('ec2')
        # Get id of current instance
        instance_id = os.popen(
            "wget -q -O - http://169.254.169.254/latest/meta-data/instance-id"
            ).read()
        # Get info of current instance from id
        instance_information = ec2.describe_instances(
            InstanceIds=[instance_id])
        # Get corresponding spot instance id
        request_id = instance_information['Reservations'][0]['Instances'][0][
            'SpotInstanceRequestId']

        return instance_id, request_id
    else:
        return None, None


def terminate_instance_and_spot_request(instance_id, spot_request_id):
    ec2 = boto3.client('ec2')

    with open("../config.json") as f:
        config = json.load(f)
        f.close()

    response = ec2.cancel_spot_instance_requests(
        DryRun=False,
        SpotInstanceRequestIds=[spot_request_id]
    )

    if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
        write_message_to_queue(
            config["aws"]["exceptions_queue"],
            "instance",
            str(response),
            None
        )
        return

    response = ec2.terminate_instances(
        DryRun=False,
        InstanceIds=[instance_id]
    )

    if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
        write_message_to_queue(
            config["aws"]["exceptions_queue"],
            "instance",
            str(response),
            None
        )
