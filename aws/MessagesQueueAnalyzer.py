import json
import pandas as pd


def summarize_execution(summary_file_path, max_epochs, out_file):
    """
    Produces an html file summarizing the state of each experiment running
    or which has run.

    Parameters
    ----------
    summary_file_path: string
        Path to the messages_queue file
    max_epochs: int
        Number of epochs the simulation is expected to run for
    out_file: string
        Path to where the output file should be saved
    """
    messages = pd.read_csv(summary_file_path, parse_dates=["timestamp"])
    experiments = messages.groupby("experiment_name")

    experiment_summaries = []

    for experiment_name in experiments:
        experiment = experiments.get_group(
            experiment_name[0]
        ).sort_values(by="timestamp", ascending=False)

        group_head = experiment.iloc[0, :]
        if len(experiment) > 1:
            time_diff = group_head["timestamp"] - experiment.iloc[1, :][
                "timestamp"]

            components = time_diff.components
            time_diff_str = "{0}h{1}m{2}s".format(
                components.hours,
                components.minutes,
                components.seconds
            )

            group_head["time_diff"] = time_diff_str
        else:
            group_head["time_diff"] = "N/A"

        percentage_progress = round(
            100 * (group_head["epoch"] * 1.) / (max_epochs - 1),
            2
        )
        group_head["progress"] = "{0}%".format(percentage_progress)
        group_head["completed"] = percentage_progress == 100 or group_head[
            "body"] in ["COMPLETED", "EXCEPTION"]

        experiment_summaries.append(group_head)

    experiment_summaries_df = pd.DataFrame(
        experiment_summaries
    ).sort_values(by="epoch", ascending=False)

    pd.set_option('display.max_colwidth', -1)
    pd.set_option('display.max_columns', None)
    with open(out_file, "w") as f:
        experiment_summaries_df.to_html(f, index=False)
        f.close()


if __name__ == "__main__":
    with open("../config.json") as f:
        config = json.load(f)
        f.close()

    summarize_execution("cloud_logs/messages_queue.csv",
                        config["num_epochs"],
                        "cloud_logs/execution_summary.html"
                        )
