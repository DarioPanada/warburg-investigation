from panaxea.core.Steppables import Helper

from aws.MessageWriter import write_message_to_queue


class HeartbeatHelper(Helper):
    """
    Simple helper that every so many epochs sends a "heartbeat" to a queue
    updating on the progress of the simulation.
    """

    def step_prologue(self, model):
        aws_config = model.properties["config"]["aws"]
        heartbeat_interval = aws_config["hearbeat_interval"]

        if model.current_epoch % heartbeat_interval == 0:
            write_message_to_queue(
                aws_config["messages_queue"],
                model.properties["name"],
                "Heartbeat at epoch {0}".format(str(model.current_epoch)),
                model
            )
