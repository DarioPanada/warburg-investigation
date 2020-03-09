import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

plt.switch_backend("agg")
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial

matplotlib.use("Qt4Agg")


def display_expected_vs_actual_curve(rep, out_path):
    """
    Given a model object. saves a scatter of the actual vs expected growth
    curve. This assumes the model has been run, otherwise no curve would be
    produced.
    Parameters
    ----------
    rep : Model
        A model instance from which growth curve values will be extracted.
    out_path : string
        The path where the image will be saved
    """
    cancer_cell_num = rep.output["agentNums"]["cancerCells"]

    cancer_cell_num = [x * 120. * 5000 / 10 ** 9 for x in cancer_cell_num]
    # plt.scatter([x/24.*2 for x in range(len(cancerCellNum))],
    # cancerCellNum, label="In-Silico")

    # Empirical data (liu et al)
    xs = [0, 13, 24]
    ys = [0, 50, 150]
    p = Polynomial.fit(xs, ys, 2)

    if len(cancer_cell_num) < 300:
        cancer_cell_num = cancer_cell_num + [cancer_cell_num[-1] for _ in
                                             range(300 - len(cancer_cell_num))]

    model_values = []

    for i, t in enumerate([x / 24. * 2 for x in range(len(cancer_cell_num))]):
        modelValue = cancer_cell_num[i]
        model_values.append(modelValue)

    cancer_cell_num = [x * 120. * 5000 / 10 ** 9 for x in cancer_cell_num]
    plt.scatter([x / 24. * 2 for x in range(len(model_values))], model_values,
                label="In-Silico")
    plt.plot(*p.linspace(), label="Empirical (Liu et al.)", color="orange")

    plt.xlabel("Time (days)")
    plt.ylabel("Volume (mm3)")
    plt.legend()
    plt.savefig(out_path)


def get_me(rep):
    """
    Given a model object, returns the mean error of the growth curve.

    Parameters
    ----------
    rep : Model
        The model instance

    Returns
    -------
    float
        The mean error on the curve
    """
    cancer_cell_num = rep.output["agentNums"]["cancerCells"]
    cancer_cell_num = [x * 120. * 5000 / 10 ** 9 for x in cancer_cell_num]
    # plt.scatter([x/24.*2 for x in range(len(cancerCellNum))],
    # cancerCellNum, label="In-Silico")

    # Empirical data (liu et al)
    xs = [0, 13, 24]
    ys = [0, 50, 150]
    p = Polynomial.fit(xs, ys, 2)

    diffs = []

    if len(cancer_cell_num) < 300:
        print("%s has %s epochs recorded, padding" % (
            rep.properties["name"], str(len(cancer_cell_num))))
        cancer_cell_num = cancer_cell_num + [cancer_cell_num[-1] for _ in
                                             range(300 - len(cancer_cell_num))]
        print(len(cancer_cell_num))

    for i, t in enumerate([x / 24. * 2 for x in range(len(cancer_cell_num))]):
        model_value = cancer_cell_num[i]
        empirical_value = p(t)
        diffs.append(abs(model_value - empirical_value))

    score = np.mean(diffs)

    return score


def get_error_series(rep, out_path):
    """
    Given a list of reports, produces a csv where each row is a simulation and
    each column the error at a certain time-point.

    Parameters
    ----------
    rep : Model
        The model instance
    out_path : string
        The path where the csv should be saved
    """
    all_errs = []
    print("hi")
    pad_count = 0
    for exp_name, k in rep.iteritems():
        cancer_cell_num = k.output["agentNums"]["cancerCells"]
        p_warburg_switch = k.properties['agents']['cancerCells'][
            'pWarburgSwitch']
        minimum_oxygen_concentration = k.properties['agents']['cancerCells'][
            'minimumOxygenConcentration']

        cancer_cell_num = [x * 120. * 5000 / 10 ** 9 for x in cancer_cell_num]
        # plt.scatter([x/24.*2 for x in range(len(cancerCellNum))],
        # cancerCellNum, label="In-Silico")

        # Empirical data (liu et al)
        xs = [0, 13, 24]
        ys = [0, 50, 150]
        p = Polynomial.fit(xs, ys, 2)

        if len(cancer_cell_num) < 300:
            cancer_cell_num = cancer_cell_num + [cancer_cell_num[-1] for _ in
                                                 range(300 - len(
                                                     cancer_cell_num))]
            print(exp_name, p_warburg_switch, minimum_oxygen_concentration)
            pad_count += 1

        errs = [exp_name, p_warburg_switch, minimum_oxygen_concentration]

        for i, t in enumerate(
                [x / 24. * 2 for x in range(len(cancer_cell_num))]):
            model_value = cancer_cell_num[i]
            empirical_value = p(t)
            err = model_value - empirical_value
            errs.append(err)

        all_errs.append(errs)
    print("hi")
    try:
        print("HERE")
        df = pd.DataFrame(all_errs)
        df.to_csv(out_path, index=False,
                  header=["name", "pWarburgSwitch",
                          "minimumOxygenConcentration"] + range(0, 300))
        print("Padded: %s" % str(pad_count))
    except Exception as e:
        print(str(e))
