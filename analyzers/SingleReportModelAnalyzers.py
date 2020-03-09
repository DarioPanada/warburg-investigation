import matplotlib.pyplot as plt
import numpy as np
import os
from panaxea.toolkit.Toolkit import depickle_from_lite

from model.agents.CancerCell import CancerCell


def get_avg_num_agents(models):
    """
    Given a set of model lite objects, returns a dictionary containing
    average number of agents at each epoch.

    Parameters
    ----------
    models : list
        The set of models lite (eg: As obtained through getModels)

    Returns
    -------
    dict
        A dictionary of average numbers of agents, with keys: cancerCells,
        tipCells, aliveCancerCells, deadCancerCells
    """
    num_cancer = np.mean(
        [m.output["agentNums"]["cancerCells"] for m in models], axis=0)
    tip_cells = np.mean([m.output["agentNums"]["tipCells"] for m in models],
                        axis=0)
    alive_cancer_cells = np.mean([m.output["agentNums"]["aliveCancerCells"]
                                  for m in models], axis=0)
    dead_cancer_cells = np.mean([m.output["agentNums"]["deadCancerCells"]
                                 for m in models], axis=0)

    return {
        "cancerCells": num_cancer,
        "tipCells": tip_cells,
        "aliveCancerCells": alive_cancer_cells,
        "deadCancerCells": dead_cancer_cells
    }


def get_avg_oxygen_concentrations(models):
    """
    Given a set of model lites, returns a dictionary containing average,
    maximum and minimum oxygen concentrations at each
    epoch.
    Parameters
    ----------
    models : list
        The set of models lite (eg: As obtained through getModels)

    Returns
    -------
    dict
        A dictionary of average, maximum and minimum oxygen concentrations.

    """
    avg_oxygen = np.mean(
        [m.output["cancerCellProperties"]["avgOxygen"] for m in models],
        axis=0)
    min_oxygen = np.mean(
        [m.output["cancerCellProperties"]["minOxygen"] for m in models],
        axis=0)
    max_oxygen = np.mean(
        [m.output["cancerCellProperties"]["maxOxygen"] for m in models],
        axis=0)

    return {
        "avg_oxygen": avg_oxygen,
        "min_oxygen": min_oxygen,
        "max_oxygen": max_oxygen
    }


def get_avg_cancer_props(models):
    """
    Given a set of model lites, returns a dictionary containing average
    cancer cell properties at each epoch.
    Parameters
    ----------
    models : list

    Returns
    -------
    dict
        The set of models lite (eg: As obtained through getModels)

    """
    avg_hif = np.mean(
        [m.output["cancerCellProperties"]["avgHif"] for m in models], axis=0)
    avg_vegf = np.mean(
        [m.output["cancerCellProperties"]["avgVegf"] for m in models], axis=0)
    avg_metabolic_rate = np.mean(
        [m.output["cancerCellProperties"]["avgMetabolicRates"] for m in
         models], axis=0)
    avg_p_synthesis = np.mean(
        [m.output["cancerCellProperties"]["avgPSynthesis"] for m in models],
        axis=0)

    return {
        "avgHif": avg_hif,
        "avgVegf": avg_vegf,
        "avgMetabolicRate": avg_metabolic_rate,
        "avgPSynthesis": avg_p_synthesis
    }


def get_avg_vegf_stimulus(models):
    """
    Given a set of model lites, returns a list of average vegf stimulus
    values at endothelial cells for each epoch.


    Parameters
    ----------
    models : list

    Returns
    -------
    list
        The aforementioned list
    """

    return np.mean(
        [m.output["endothelialCellProperties"]["avgVegf"] for m in models],
        axis=0)


def get_avg_tumour_volume(models):
    """
    Given a set of model lites, returns a list of average tumour volumes for
    each epoch.

    Parameters
    ----------
    models : list
        The set of models lite (eg: As obtained through getModels)

    Returns
    -------
    list
        The aforementioned list
    """

    return np.mean([m.output["maxDistances"] for m in models], axis=0)


def get_avg_num_warburg_cells(models):
    """
    Given a set of model lites, returns a list of PERCENTAGE fo warburg
    cells for each epoch.

    Parameters
    ----------
    models : list
        The set of models lite (eg: As obtained through getModels)

    Returns
    -------
    list
        The aforementioned list
    """
    return np.mean(
        [m.output["cancerCellProperties"]["numWarburgCells"] for m in models],
        axis=0)


def get_post_execution_analysis(target_dir, imgs_dir="imgs"):
    """
    Given a set of pickle lite objects, one for each end-state of models (eg
    as generated via the runModels function) this
    runs all post-execution analysis to generate relevant scatters and saves
    them to a specified directory. Average
    values across all execution are used to generate these EXCEPT FOR (see
    below):

    Parameters
    ----------
    target_dir : string
        The report dir, this is returned by the runModels function
    imgs_dir : string
        The directory where the images will be saved (it is created if it
        does not exist)
    """

    imgs_path = "%s/%s" % (target_dir, imgs_dir)
    if not os.path.exists(imgs_path):
        os.mkdir(imgs_path)

    pickles = [f for f in os.listdir(target_dir) if f.endswith(".pickle")]

    if len(pickles) != 1:
        print("Expecting one pickle object, found {0}: {1}".format(
            len(pickles), ",".join(pickles)))

    model = depickle_from_lite("{0}/{1}".format(target_dir, pickles[0]))

    c = CancerCell(model)
    get_hif_from_oxygen(c, out_path=imgs_path, interval=0.001)
    get_metabolic_rate_from_hif(c, out_path=imgs_path, interval=0.001)
    get_probability_synthesis_from_hif(c, out_path=imgs_path, interval=0.001)
    get_vegf_secretion_rate_from_hif_concentration(
        c,
        out_path=imgs_path,
        interval=0.001)
    c.warburg_switch = True
    get_hif_from_oxygen(
        c,
        out_path=imgs_path,
        interval=0.001,
        suffix="warburg")

    models = [model]
    post_execution_agent_num_visualizer(
        render=False,
        num_agents=get_avg_num_agents(models),
        out_path=imgs_path)
    post_execution_cancer_cell_properties_visualizer(
        render=False,
        avg_props=get_avg_cancer_props(models),
        out_path=imgs_path)
    vegf_stimulus_viewer(
        render=False,
        avg_vegf_stimulus=get_avg_vegf_stimulus(models),
        out_path=imgs_path)
    post_execution_oxygen_concentration_visualizer(
        render=False,
        avgProps=get_avg_oxygen_concentrations(models),
        out_path=imgs_path)
    tumour_volume_viewer(
        render=False,
        avg_tumour_volume=get_avg_tumour_volume(models),
        out_path=imgs_path)
    warburg_num_viewer(
        render=False,
        avg_warburg_cells=get_avg_num_warburg_cells(models),
        out_path=imgs_path)

    visualize_glucose_distributions(
        model.output["cancerCellProperties"]["GlucoseDistributions"],
        out_path=imgs_path)
    visualize_oxygen_distributions(
        model.output["cancerCellProperties"]["OxygenDistributions"],
        out_path=imgs_path)

    visualize_final_summary_cancer_cell_death(
        model.output["causesOfDeath"],
        out_path=imgs_path)
    output_avg_age_cell_death(
        model.output["causesOfDeath"],
        out_path=imgs_path)
    save_oxygen_distributions(
        model.output["cancerCellProperties"]["OxygenDistributions"], imgs_path)
    save_glucose_distributions(
        model.output["cancerCellProperties"]["GlucoseDistributions"],
        imgs_path)
    save_hif_distributions(model.output["cancerCellProperties"][
                               "HIFExpressionRatesDistributions"], imgs_path)

    oxygen_dists = model.output["cancerCellProperties"]["OxygenDistributions"]
    hif_dists = model.output["cancerCellProperties"][
        "HIFExpressionRatesDistributions"]

    if len(oxygen_dists) > 0 and len(hif_dists) > 0:
        save_final_oxygen_hif_distributions(oxygen_dists, hif_dists, imgs_path)


def visualize_glucose_distributions(glucose_distributions, out_path):
    """
    Creates a scatter used to show average, maximum and minimum glucose
    concentrations for each epoch.

    Parameters
    ----------
    glucose_distributions : list
        An appropriate data-structure of concentrations as created by
        GlucoseConcentrationWatcher
helper
    out_path : string
        The location where the file will be saved
    """
    maxes = []
    mins = []
    avgs = []

    for i, d in enumerate(glucose_distributions):
        mn = d["bins"][0]
        mx = d["bins"][-1]

        accumulator = 0

        for i in range(len(d["n"])):
            accumulator += d["n"][i] * d["bins"][i]

        avg = accumulator / sum(d["n"])

        avgs.append(avg)

        mins.append(mn)
        maxes.append(mx)

    plt.figure()
    plt.scatter(range(len(maxes)), maxes, label="max")
    plt.scatter(range(len(maxes)), mins, label="min")
    plt.scatter(range(len(maxes)), avgs, label="avg")
    plt.xlabel("Epoch")
    plt.ylabel("Glucose")
    plt.title("Glucose in Tumour")
    plt.legend()
    plt.savefig("%s/glucose_concentrations_timeline.png" % out_path)


def visualize_oxygen_distributions(oxygen_distributions, out_path):
    """
    Creates a scatter used to show average, maximum and minimum oxygen
    concentrations for each epoch.
    Parameters
    ----------
    oxygen_distributions : object
        An appropriate data-structure of concentrations as created by
        OxygenConcentrationWatcher helper
    out_path : string
        The location where the file will be saved
    """
    maxes = []
    mins = []
    avgs = []

    for i, d in enumerate(oxygen_distributions):
        mn = d["bins"][0]
        mx = d["bins"][-1]

        accumulator = 0

        for i in range(len(d["n"])):
            accumulator += d["n"][i] * d["bins"][i]

        avg = accumulator / sum(d["n"])

        avgs.append(avg)

        mins.append(mn)
        maxes.append(mx)

    plt.figure()
    plt.scatter(range(len(maxes)), maxes, label="max")
    plt.scatter(range(len(maxes)), mins, label="min")
    plt.scatter(range(len(maxes)), avgs, label="avg")
    plt.xlabel("Epoch")
    plt.ylabel("Oxygen (mmHg)")
    plt.title("Oxygen in Tumour")
    plt.legend()
    plt.savefig("%s/oxygen_concentrations_timeline.png" % out_path)


def output_avg_age_cell_death(causes_of_death, out_path):
    """
    For each category of cell (warburgDeathGlucose, warburgDeathOxygen,
    nonWarburgDeathGlucose, nonWarburgDeathOxygen) it prints the average age
    at
    which the cell died (and standard deviation).

    Parameters
    ----------
    causes_of_death : list
        Causes of Death array, as created by DeathCauseWatcher helper
    out_path : string
        Location where the file will be saved
    """
    causes_of_death = causes_of_death[-1]

    with open("%s/ageOfDeath.csv" % out_path, "w") as f:
        f.write("class,avg,stDev\n")
        f.write("warburgGlucose,%s,%s\n" % (
            causes_of_death["warburgDeathGlucose"]["avgAge"],
            causes_of_death["warburgDeathGlucose"]["stDev"]))
        f.write("warburgOxygen,%s,%s\n" % (
            causes_of_death["warburgDeathOxygen"]["avgAge"],
            causes_of_death["warburgDeathOxygen"]["stDev"]))
        f.write("nonWarburgGlucose,%s,%s\n" % (
            causes_of_death["nonWarburgDeathGlucose"]["avgAge"],
            causes_of_death["nonWarburgDeathGlucose"]["stDev"]))
        f.write("nonwarburgOxygen,%s,%s\n" % (
            causes_of_death["nonWarburgDeathOxygen"]["avgAge"],
            causes_of_death["nonWarburgDeathOxygen"]["stDev"]))
        f.close()


def visualize_final_summary_cancer_cell_death(causes_of_death, out_path):
    """
    Visualizes, for the last epoch, the cumulative sum of cancer cell death
    cause.
    Parameters
    ----------
    causes_of_death : list
        A list of dictionaries bound to the causesOfDeath model property,
        as generated by DeathCauseWatcher helper
    out_path : string
        Directory where the pictures will be saved
    """
    causes_of_death = causes_of_death[-1]

    warburg_glucose = causes_of_death["warburgDeathGlucose"]["num"]
    warburg_oxygen = causes_of_death["warburgDeathOxygen"]["num"]
    non_warburg_oxygen = causes_of_death["nonWarburgDeathOxygen"]["num"]
    non_warburg_glucose = causes_of_death["nonWarburgDeathGlucose"]["num"]

    plt.figure()
    plt.title("Causes of Death")
    plt.ylabel("Number of Agents")
    labels = ["warburg_glucose", "non_warburg_oxygen", "non_warburg_glucose",
              "warburg_oxygen"]
    plt.bar([0, 1, 2, 3],
            [warburg_glucose, non_warburg_oxygen, non_warburg_glucose,
             warburg_oxygen])
    plt.xticks([0, 1, 2, 3], labels)
    plt.savefig("%s/causesOfDeath.png" % out_path)
    plt.close()


def save_final_oxygen_hif_distributions(oxygen_distributions,
                                        hif_distributions, imgs_path):
    """
    Saves final stacked bar charts representing end-state oxygen and hif
    distributions.

    Parameters
    ----------
    oxygen_distributions : list
        Oxygen distributions at end
    hif_distributions : list
        HIF distributions at end
    imgs_path : string
        Path to save image
    """
    fig = plt.figure()
    fig.add_subplot(1, 2, 1)

    concentrations_at_end = hif_distributions
    concentrations_zipped = [
        (concentrations_at_end["bins"][i], concentrations_at_end["n"][i]) for i
        in range(len(concentrations_at_end["bins"]) - 1)]

    basic = [c for c in concentrations_zipped if c[0] < 2]
    basic = sum(b[1] for b in basic)
    step = 2
    intervals = np.arange(2, 16, step)
    enhanceds = []
    for i in intervals:
        enhanceds.append((
            "%s-%s" % (str(round(i, 2)), str(round(i + step, 2))),
            sum(c[1] for c in concentrations_zipped if
                i <= c[0] < i + step)))
    total = basic + sum(e[1] for e in enhanceds)

    plt.bar([0], basic / total)
    enhanceds_percentage = [(e[0], e[1] / total) for e in enhanceds]
    accumulator = 0
    plt.bar([1], enhanceds_percentage[0][1], label=enhanceds[0][0])
    accumulator += enhanceds_percentage[0][1]

    del enhanceds[0]

    for e in enhanceds_percentage:
        plt.bar([1], e[1], label=e[0], bottom=accumulator)
        accumulator += e[1]

    plt.legend()
    plt.xticks([0, 1], ["base", "enhanced"])
    plt.title("Distribution of HIF Expression Rates")
    plt.ylabel("Number of Cells (Percentage)")
    totalOxygen = total

    concentrations_at_end = oxygen_distributions
    concentrations_zipped = [
        (concentrations_at_end["bins"][i], concentrations_at_end["n"][i]) for i
        in range(len(concentrations_at_end["bins"]) - 1)]
    fig.add_subplot(1, 2, 2)

    hypoxics = []
    step = 0.05
    intervals = np.arange(0, 0.2, step)
    for i in intervals:
        hypoxics.append(("%s%% to %s%%" % (
            str(round(i, 2)), str(round(i + step, 2))), sum(
            c[1] for c in concentrations_zipped if i <= c[0] < i + step)))

    normoxics = [c[1] for c in concentrations_zipped if c[0] >= i + step]

    tot_hypoxics = sum(h[1] for h in hypoxics)
    tot_normoxics = sum(normoxics)
    total = tot_normoxics + tot_hypoxics

    hypoxics = [(h[0], h[1] / total) for h in hypoxics]

    plt.bar([0], tot_normoxics / total)

    plt.bar([1], hypoxics[0][1], label=hypoxics[0][0])
    prev = hypoxics[0][1]

    del hypoxics[0]

    for h in hypoxics:
        plt.bar([1], h[1], label=h[0], bottom=prev)
        prev += h[1]

    plt.xticks([0, 1], ["normoxic", "hypoxic"])
    plt.ylabel("Number of Cells (Percentage)")
    plt.title("Distribution of Oxygen Concentrations")
    plt.legend()
    plt.text(-3, 0, "Total Agents Oxygen %s - Total Agents HIF %s" % (
        str(totalOxygen), str(total)))
    plt.savefig("%s/oxygen_hif_distributions_end.png" % imgs_path)


def save_glucose_distributions(distributions, imgs_path,
                               distributions_folder_name="glucoseDistributions"
                               ):
    """
    From a single model loads and saves figures representing glucose
    concentration distributions
    Parameters
    ----------
    distributions : string
        List of distributions generated by GlucoseConcentratrionWatcher
    imgs_path : string
        Directory where subdirectory will be created
    distributions_folder_name : string, optional
        name of directory where figures will be saved. Defaults to
        glucoseDistributions
    """
    distributions_dir = "%s/%s" % (imgs_path, distributions_folder_name)

    if not os.path.exists(distributions_dir):
        os.mkdir(distributions_dir)

    for d in distributions:
        plt.figure()
        n_percentage = [n / sum(d["n"]) for n in d["n"]]

        plt.bar(d["bins"][:-1], n_percentage, width=1)
        plt.xlabel("GlucoseConcentration")
        plt.ylabel("Amount of Cancer Cells(Percentage)")
        plt.title("Glucose Concentration in Tumour")
        plt.savefig("%s/distribution_%s.png" % (distributions_dir, d["epoch"]))
        plt.close()


def save_oxygen_distributions(distributions, imgs_path,
                              distributions_folder_name="oxygenDistributions"):
    """
    From a single model loads and saves figures representing oxygen
    concentration distributions

    Parameters
    ----------
    distributions : list
        List of distributions generated by OxygenConcentratrionWatcher
    imgs_path : string
        Directory where subdirectory will be created
    distributions_folder_name : string, optional
        name of directory where figures will be saved. Defaults to
        oxygenDistributions
    """
    distributions_dir = "%s/%s" % (imgs_path, distributions_folder_name)

    if not os.path.exists(distributions_dir):
        os.mkdir(distributions_dir)

    for d in distributions:
        plt.figure()
        n_percentage = [n / sum(d["n"]) for n in d["n"]]

        plt.bar(d["bins"][:-1], n_percentage, width=1)
        plt.xlim([d["bins"][0] - 1, d["bins"][-1] + 1])
        plt.xlabel("PPO2 (mmHg)")
        plt.ylabel("Amount of Cancer Cells(Percentage)")
        plt.title("Oxygen Concentration in Tumour")
        plt.savefig("%s/distribution_%s.png" % (distributions_dir, d["epoch"]))
        plt.close()


def save_hif_distributions(distributions, imgs_path,
                           distributions_folder_name="HIFDistributions"):
    """
    From a single model loads and saves figures representing hif
    concentration distributions

    Parameters
    ----------
    distributions : list
        List of distributions generated by CancerCellWatcher
    imgs_path : string
        Directory where subdirectory will be created
    distributions_folder_name : string, optional
        Name of directory where figures will be saved. Optional, defaults to
        HIFDistributions
    """
    distributions_dir = "%s/%s" % (imgs_path, distributions_folder_name)

    if not os.path.exists(distributions_dir):
        os.mkdir(distributions_dir)

    for d in distributions:
        plt.figure()
        n_percentage = [n / sum(d["n"]) for n in d["n"]]

        plt.bar(d["bins"][:-1], n_percentage, width=1)
        plt.xticks(np.arange(1, 16, 1))
        plt.xlabel("HIF Expression Rate")
        plt.ylabel("Amount of Cancer Cells(Percentage)")
        plt.title("HIF Expression Rates in Tumour")
        plt.savefig("%s/distribution_%s.png" % (distributions_dir, d["epoch"]))
        plt.close()


def post_execution_agent_num_visualizer(model=None, num_agents=None,
                                        render=True, out_path=None):
    """
    Creates a scatter plot summarizing the amount of each agent type at each
    epoch.

    Parameters
    ----------
    model : Model
        The model object
    num_agents : A dictionary containing agent numbers at each epoch
    render : bool, optional
        Defaults to true, if true the scatters are displayed
    out_path : string, optional
        Defaults ot none, if not none the figure is saved to
        out_path/numAgents.png, or a path can be specified
    """
    if model is None and num_agents is None:
        print("Need to provide model or dictionary of agents nums")
        exit()

    if model is not None:
        num_cancer = model.output["agentNums"]["cancerCells"]
        tip_cells = model.output["agentNums"]["tipCells"]
        alive_cancer_cells = model.output["agentNums"]["aliveCancerCells"]
        dead_cancer_cells = model.output["agentNums"]["deadCancerCells"]
    else:
        num_cancer = num_agents["cancerCells"]
        tip_cells = num_agents["tipCells"]
        alive_cancer_cells = num_agents["aliveCancerCells"]
        dead_cancer_cells = num_agents["deadCancerCells"]

    x_epochs = range(len(num_cancer))

    plt.figure()
    plt.scatter(x_epochs, num_cancer, label="Total Cancer Cells")
    plt.scatter(x_epochs, alive_cancer_cells, label="Alive Cancer Cells")
    plt.scatter(x_epochs, dead_cancer_cells, label="Dead Cancer Cells")
    plt.legend()

    plt.xlabel("Epoch")
    plt.ylabel("Number of Cancer Cells")
    plt.title("Cancer Size in Time")

    if out_path is not None:
        plt.savefig("%s/numCancerCells.png" % out_path)

    plt.figure()
    plt.scatter(x_epochs, tip_cells)
    plt.xlabel("Epoch")
    plt.ylabel("Number of Endothelial Cells")
    plt.title("Blood Vessel Development in Time")

    if out_path is not None:
        plt.savefig("%s/numEndothelials.png" % out_path)

    if render:
        plt.show()


def post_execution_oxygen_concentration_visualizer(model=None, avgProps=None,
                                                   render=True, out_path=None):
    """
    Creates a scatter for average, max and min oxygen concentrations in time.

    Parameters
    ----------
    model : Model
        The model object
    avgProps : dict
        Average oxygen properties
    render : bool, optional
        Defaults to true, if true the scatters are displayed
    out_path : string, optional
         Defaults ot none, if not none the figure is saved to
         out_path/numAgents.png
    """

    if model is not None:
        avg_oxygen = model.output["cancerCellProperties"]["avgOxygen"]
        max_oxygen = model.output["cancerCellProperties"]["maxOxygen"]
        min_oxygen = model.output["cancerCellProperties"]["minOxygen"]
    elif avgProps is not None:
        avg_oxygen = avgProps["avg_oxygen"]
        min_oxygen = avgProps["min_oxygen"]
        max_oxygen = avgProps["max_oxygen"]
    else:
        print("Need to provide model or dictionary of oxygen properties")
        exit()

    x_epochs = range(len(avg_oxygen))

    fig = plt.figure()

    fig.add_subplot(1, 3, 1)
    plt.scatter(x_epochs, avg_oxygen)
    plt.title("Average Oxygen in Time")
    plt.xlabel("Epoch")
    plt.ylabel("Oxygen Concentration")

    fig.add_subplot(1, 3, 2)
    plt.scatter(x_epochs, max_oxygen)
    plt.title("Max Oxygen in Time")
    plt.xlabel("Epoch")
    plt.ylabel("Oxygen Concentration")

    fig.add_subplot(1, 3, 3)
    plt.scatter(x_epochs, min_oxygen)
    plt.title("Minimum Oxygen in Time")
    plt.xlabel("Epoch")
    plt.ylabel("Oxygen Concentration")

    if out_path is not None:
        fig.savefig("%s/oxygenConcentrations.png" % out_path)

    if render:
        plt.show()


def post_execution_cancer_cell_properties_visualizer(model=None,
                                                     avg_props=None,
                                                     render=True,
                                                     out_path=None):
    """
    Creates a scatter showing the properteis of cancer cells in time
    Parameters
    ----------
    model : Model
        The model object
    avg_props : dict
        A dictionary containing average cancer cell properties at each epoch
    render : bool, optional
        Defaults to true, if true the scatters are displayed
    out_path : string, optional
        Defaults ot none, if not none the figure is saved to
        out_path/numAgents.png
    """
    if model is None and avg_props is None:
        print("Need to provide model or dictionary of agents nums")
        exit()

    if model is not None:
        avg_hif = model.output["cancerCellProperties"]["avgHif"]
        avg_vegf = model.output["cancerCellProperties"]["avgVegf"]
        avg_metabolic_rate = model.output["cancerCellProperties"][
            "avgMetabolicRate"]
        avg_p_synthesis = model.output["cancerCellProperties"]["avgPSynthesis"]
    else:
        avg_hif = avg_props["avgHif"]
        avg_vegf = avg_props["avgVegf"]
        avg_metabolic_rate = avg_props["avgMetabolicRate"]
        avg_p_synthesis = avg_props["avgPSynthesis"]

    x_epochs = range(len(avg_hif))

    rows = cols = 2
    fig = plt.figure()

    x_template = "Avg %s"
    y_template = "Avg %s in time"

    fig.add_subplot(rows, cols, 1)
    plt.scatter(x_epochs, avg_hif)
    plt.xlabel("Epoch")
    target = "HIF Expression"
    plt.ylabel(x_template % target)
    plt.title(y_template % target)

    fig.add_subplot(rows, cols, 2)
    plt.scatter(x_epochs, avg_vegf)
    plt.xlabel("Epoch")
    target = "VEGF Expression Rate"
    plt.ylabel(x_template % target)
    plt.title(y_template % target)

    fig.add_subplot(rows, cols, 3)
    plt.scatter(x_epochs, avg_metabolic_rate)
    plt.xlabel("Epoch")
    target = "Metabolic Rate"
    plt.ylabel(x_template % target)
    plt.title(y_template % target)

    fig.add_subplot(rows, cols, 4)
    plt.scatter(x_epochs, avg_p_synthesis)
    plt.xlabel("Epoch")
    target = "Probability Progress in Synthesis"
    plt.ylabel(x_template % target)
    plt.title(y_template % target)

    if out_path is not None:
        fig.savefig("%s/cancerCellProperties.png" % out_path)

    if render:
        plt.show()


def get_hif_from_oxygen(cancer_cell, interval=10, render=False, out_path=None,
                        suffix=""):
    """
    Given a cancer cell object, calculates hif expression rates for all
    oxygen PPO2 between 0 and 200

    Parameters
    ----------
    cancer_cell : CancerCell
        The cancer cell object
    interval : float, optional
        The x-interval (eg setting this to 0.1 will calculate for 0, 0.1,
        0.2, ..., 0.9, 1) Defaults to 10
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    out_path : bool, optional
         If set to true, saves the figure to this path
    suffix : string, optional
        If set, will add a suffix to the end of the file name. Defaults to the
        empty string

    Returns
    -------
    list
        An array of tuples (oxygenConcentratoin, hif expression)

    """
    concentrations = np.arange(0, 200, interval)
    if not cancer_cell.warburg_switch:
        hif_expressions = [
            cancer_cell._calculate_hif_expression_rate_from_oxygen(c)
            for c in concentrations]
    else:
        hif_expressions = [
            cancer_cell._calculate_hif_expression_rate_from_oxygen_warburg(c)
            for c in concentrations]

    plt.figure()
    plt.scatter(concentrations, hif_expressions)
    plt.xlabel("Oxygen Concentration")
    plt.ylabel("Hif Expression Rate")
    plt.title("HIF Expression Rate Across Multiple oxygen Concentrations")

    if out_path is not None:
        plt.savefig("%s/oxygen_to_hif_%s.png" % (out_path, suffix))

    if render:
        plt.show()
    return zip(concentrations, hif_expressions)


def get_metabolic_rate_from_hif(cancer_cell, interval=1, render=False,
                                out_path=None):
    """
    Given a cancer cell object, calculates metabolic rate for all hif
    expression rates between 1 and 15

    Parameters
    ----------
    cancer_cell : CancerCell
        The cancer cell object
    interval : float, optional
        The x-interval (eg setting this to 1 will calculate for 0, 1,2,3,
        ...,14, 15) Defaults to 1
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    out_path : string, optional
         If set to a value, saves the figure to this path. Defaults to None
    """
    hif_rates = np.arange(1, 15, interval)

    def get_metabolic_rate(h, a):
        a.currentHifRate = h
        a._update_metabolic_rate()
        mr = a.current_metabolic_rate

        return mr

    metabolic_rates = [get_metabolic_rate(h, cancer_cell) for h in hif_rates]

    plt.figure()
    plt.scatter(hif_rates, metabolic_rates)
    plt.xlabel("HIF Expression Rate")
    plt.ylabel("Metabolic Rate")
    plt.legend()
    plt.title("Change in Metabolic Rate across HIF Expression Rates")

    if out_path is not None:
        plt.savefig("%s/hif_to_metabolic_rate.png" % out_path)

    if render:
        plt.show()

    return zip(hif_rates, metabolic_rates)


def get_probability_synthesis_from_hif(cancer_cell, interval=1, render=False,
                                       out_path=None):
    """
    Given a cancer cell object, calculates the probability of the agent
    progressing into synthesis for all hif expression
    rates between 1 and 15

    Parameters
    ----------
     cancer_cell : CancerCell
        The cancer cell object
    interval : float, optional
        The x-interval (eg setting this to 1 will calculate for 0, 1,2,3,
        ...,14, 15) Defaults to 1
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    out_path : string, optional
         If set to a value, saves the figure to this path. Defaults to None
    """
    hif_rates = np.arange(1, 15, interval)

    def get_p_synthesis(h, a):
        a.currentHifRate = h
        a._update_p_synthesis()
        mr = a.current_p_synthesis

        return mr

    probabilities = [get_p_synthesis(h, cancer_cell) for h in hif_rates]

    plt.figure()
    plt.scatter(hif_rates, probabilities)
    plt.xlabel("HIF Expression Rate")
    plt.ylabel("Probability Synthesis")
    plt.legend()
    plt.title(
        "Change in Probability of Progressing into Synthesis across HIF "
        "Expression Rates")

    if out_path is not None:
        plt.savefig("%s/hif_to_p_synthesis.png" % out_path)

    if render:
        plt.show()

    return zip(hif_rates, probabilities)


def get_vegf_secretion_rate_from_hif_concentration(cancer_cell, interval=1,
                                                   render=False,
                                                   out_path=None):
    """
    Given a cancer cell object, calculates its vegf secretion rate for all
    hif expression rates between 1 and 15

    Parameters
    ----------
     cancer_cell : CancerCell
        The cancer cell object
    interval : float, optional
        The x-interval (eg setting this to 1 will calculate for 0, 1,2,3,
        ...,14, 15) Defaults to 1
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    out_path : string, optional
         If set to a value, saves the figure to this path. Defaults to None
    """
    hif_rates = np.arange(1, 15, interval)

    def get_vegf_secretion_rate(h, a):
        a.currentHifRate = h
        a._update_vegf_secretion_rate()
        mr = a.current_vegf_secretion_rate

        return mr

    rates = [get_vegf_secretion_rate(h, cancer_cell) for h in hif_rates]

    plt.figure()
    plt.scatter(hif_rates, rates)
    plt.xlabel("HIF Expression Rate")
    plt.ylabel("VEGF Expression Rate")
    plt.legend()
    plt.title("Change in VEGF ExpressionRates across HIF Expression Rates")

    if out_path is not None:
        plt.savefig("%s/hig_to_vegf_secretion_rate.png" % out_path)

    if render:
        plt.show()

    return zip(hif_rates, rates)


def warburg_num_viewer(avg_warburg_cells, out_path, render=False):
    """
    Displays a scatter showing the average number of warburg cells in time
    Parameters
    ----------
    avg_warburg_cells : list
        Average number of warburg cells at each epoch.
    out_path : string, optional
        If set to a value, saves the figure to this path. Defaults to None
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    """
    x_epochs = range(len(avg_warburg_cells))

    plt.figure()
    plt.scatter(x_epochs, avg_warburg_cells)
    plt.xlabel("Epoch")
    plt.ylabel("Avg Number of Warburg Cells (Percentage)")
    plt.title("Average Number of Warburg Cells in Time")

    if out_path is not None:
        plt.savefig("%s/warburgCells.png" % out_path)

    if render:
        plt.show()


def tumour_volume_viewer(model=None, avg_tumour_volume=None, render=False,
                         out_path=None):
    """
    Produces a scatter showing tumour volume in time.

    Parameters
    ----------
    model : Model
        The model object
    avg_tumour_volume : dict
        A list containing average tumour volume
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    out_path : string, optional
         If set to a value, saves the figure to this path. Defaults to None
    """
    if model is not None:
        avg_tumour_volume = model.output["maxDistances"]
    x_epochs = range(len(avg_tumour_volume))

    plt.figure()
    plt.scatter(x_epochs, avg_tumour_volume)
    plt.xlabel("Epoch")
    plt.ylabel("Avg Tumour Volume")
    plt.title("Average Tumour Volume in Time")

    if out_path is not None:
        plt.savefig("%s/tumourVolume.png" % out_path)

    if render:
        plt.show()


def vegf_stimulus_viewer(model=None, avg_vegf_stimulus=None, render=True,
                         out_path=None):
    """
    Produces a scatter showing average VEGF stimulus in time.

    Parameters
    ----------
    model : Model
        The model object
    avg_tumour_vegf_stimulus : dict
        A list containing average vegf stimulus in time
    render : bool, optional
        If set to true, a scatter of values is displayed, defaults to false
    out_path : string, optional
         If set to a value, saves the figure to this path. Defaults to None
    """

    if model is None and avg_vegf_stimulus is None:
        print("You must provide model or vegfStimulus")
        exit()

    if model is not None:
        avg_vegf_stimulus = model.output["endothelialCellProperties"][
            "avgVegf"]

    x_epochs = range(len(avg_vegf_stimulus))

    plt.figure()
    plt.scatter(x_epochs, avg_vegf_stimulus)
    plt.xlabel("Epoch")
    plt.ylabel("Avg VEGF Concentration")
    plt.title("Average VEGF Concentration at Blood Vessels in Time")

    if out_path is not None:
        plt.savefig("%s/vegfStimulus.png" % out_path)

    if render:
        plt.show()
