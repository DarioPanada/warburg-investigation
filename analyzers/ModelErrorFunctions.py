import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
matplotlib.use("Agg")

plt.switch_backend("agg")
from numpy.polynomial import Polynomial
from panaxea.toolkit.Toolkit import depickle_from_lite

matplotlib.use("Qt4Agg")

"""
This script contains functions to calculate and visualize expected vs actual
growth curves and calculate errors.
"""


def convert_num_agents_to_volume(num_cancer_cell_agents,
                                 cancer_cells_per_agent,
                                 cancer_cell_volume):
    """
    Given an amount of cancer cell agents, converts this to the
    corresponding volume.

    Parameters
    ----------
    num_cancer_cell_agents : int
        The number of cancer cell agents present at a given time
    cancer_cells_per_agent : int
        Number of cancer cells we assume each agent represents
    cancer_cell_volume : int
        Volume of a cancer cell - This should be provided in um3

    Returns
    -------
    float
        The calculated volume - in mm3
    """

    num_cancer_cells = num_cancer_cell_agents * cancer_cells_per_agent
    volume_mm_3 = num_cancer_cells * cancer_cell_volume
    volume_um_3 = (1. * volume_mm_3) / 10 ** 9
    return volume_um_3


def get_fitness_function_polynomial():
    """
    Returns a polynomial function p that can be used to calculate the expected
    tumour volume at a time t, 0 <= t <= 300. The time is expressed in DAYS.

    Eg: p(2) returns the expected volume at day 2.

    Returns
    -------
    Polynomial
        The polynomial function.
    """
    xs = [0, 13, 24]
    ys = [0, 50, 150]
    p = Polynomial.fit(xs, ys, 2)

    return p


def epoch_to_day(epoch, epoch_duration):
    """
    Given an epoch, calculates the real-time duration in days. (This may
    return a floating point value, such as day 14.73)
    Parameters
    ----------
    epoch : int
        The epoch we wish to convert
    epoch_duration : int
        The duration of a single epoch, in hours

    Returns
    -------
    float
        The corresponding day value
    """

    return round(epoch * epoch_duration / 24., 2)


def get_cancer_volume_series_from_agent_nums(cancer_cell_num, max_epochs):
    """
    From a series of cancer cell numbers returns a series of calculated
    tumour volumes in time
    Parameters
    ----------
    cancer_cell_num : list
        List of integers reporting the number of cancer cell agents in time
    max_epochs : int
        Maximum number of epochs the simulation is expected to run for
    Returns
    -------
    list
        List where each element reports the calculated volume at such an epoch
    """
    num_epochs = len(cancer_cell_num)

    if num_epochs < max_epochs:
        missing_epochs = max_epochs - num_epochs
        padding = [cancer_cell_num[-1]] * missing_epochs
        cancer_cell_num = cancer_cell_num + padding

    cancer_volumes = [convert_num_agents_to_volume(n, 120, 5000)
                      for n in cancer_cell_num]

    return cancer_volumes


def get_expected_volume_series(max_epochs, epoch_duration):
    """
    Returns a list where each element corresponds to the expected volume.

    Values are returned in m33

    Parameters
    ----------
    max_epochs : int
        Maximum number of epochs the simulation is expected to run for
    epoch_duration : int
        Number of hours per epoch

    Returns
    -------
    list
        List where each element reports the expected volume at such an
        epoch, in mm3
    """

    p = get_fitness_function_polynomial()
    day_values = [epoch_to_day(e, epoch_duration) for e in range(max_epochs)]
    expected_volumes = [p(d) for d in day_values]

    return expected_volumes


def get_error_series(rep, max_epochs, epoch_duration):
    """
    Given a model object, returns the error at each epoch intended as the
    difference between expected and actual growth curves.

    So, a positive error means the model reported a smaller volume than
    expected, whereas a negative error means the model reported a larger
    volume than expected.

    Parameters
    ----------
    rep : Model
        The model instance
    epoch_duration : int
        Number of hours per epoch
    max_epochs : int
        The maximum number of epochs the model is expected to run for

    Returns
    -------
    float
        The mean error on the curve
    """
    # Number of cancer cells at each epoch
    cancer_cell_num = rep.output["agentNums"]["cancerCells"]

    cancer_volumes = get_cancer_volume_series_from_agent_nums(cancer_cell_num,
                                                              max_epochs)

    expected_volumes = get_expected_volume_series(max_epochs, epoch_duration)

    errors = [expected_volumes[e] - cancer_volumes[e]
              for e in range(max_epochs)]

    return errors


def visualize_actual_vs_expected_curves(rep, max_epochs, epoch_duration,
                                        show_figure=True, out_path=None):
    """
    Creates a scatter displaying actual vs expected growth curves for a
    model report.
    Parameters
    ----------
    rep : Model
        The model instance
    max_epochs : int
        The maximum number of epochs the simulation is supposed to run for
    epoch_duration : int
        Number of hours per epoch
    show_figure : bool, optional
        If set to true, the scatter is displayed. Defaults to true.
    out_path : string, optional
        If not None, this should be a full path to a file where the scatter
        should be saved. Defaults to None
    """
    cancer_cell_num = rep.output["agentNums"]["cancerCells"]

    cancer_volumes = get_cancer_volume_series_from_agent_nums(cancer_cell_num,
                                                              max_epochs)

    expected_volumes = get_expected_volume_series(max_epochs, epoch_duration)

    plt.figure()

    xs = range(max_epochs)

    plt.scatter(xs, cancer_volumes, label="Model Growth")
    plt.scatter(xs, expected_volumes, label="Expected Growth")

    plt.title("Model vs Expected Growth")

    day_labels = [str(epoch_to_day(e, epoch_duration))
                  if epoch_to_day(e, epoch_duration) % 5 == 0
                  else "" for e in xs]
    plt.xticks(xs, day_labels)

    plt.xlabel("Time (days)")
    plt.ticklabel_format(axis="y", style="sci")
    plt.ylabel("Volume (mm3)")
    plt.legend()

    if show_figure:
        plt.show()

    if out_path is not None:
        plt.savefig(out_path)

    plt.close()


if __name__ == "__main__":
    # Sample error calculation and displaying

    max_epochs = 300
    epoch_duration = 2

    reports_dir = "../reports/"
    experiment_dirs = [d for d in os.listdir(reports_dir)
                       if os.path.isdir(reports_dir + d)]
    print("The following {0} experiment directories were found: {1}".format(
        len(experiment_dirs),
        experiment_dirs))

    # Going through each experiment
    for experiment_dir in experiment_dirs:
        full_path = reports_dir + experiment_dir

        # Assuming one pickle file per directory
        pickle = [reports_dir + experiment_dir + "/" + f
                  for f in os.listdir(full_path)
                  if f.endswith(".pickle")][0]
        model = depickle_from_lite(pickle)

        visualize_actual_vs_expected_curves(
            model,
            max_epochs,
            epoch_duration,
            show_figure=False,
            out_path=experiment_dir + ".jpg")

        error_series = get_error_series(model, max_epochs, epoch_duration)
        abs_errors = map(abs, error_series)
        mean_error = np.mean(abs_errors)
        print(experiment_dir)
        print("Mean Error: {0} mm3".format(round(mean_error, 2)))
