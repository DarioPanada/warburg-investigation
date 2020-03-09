import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.switch_backend("agg")
import numpy as np
from panaxea.core.Steppables import Helper


class OxygenConcentrationWatcher(Helper, object):
    """
    Keeps track of average, minimum and maximum oxygen concentration at cancer
    cells at each epoch.

    Output keys are created under "cancerCellProperties".

    Attributes
    ----------
    model : Model
        The model instance
     cancer_cell_class_name : string, optional
        The name of the cancer cell class, defaults to Cancer Cell
    interval : int
        The number of epochs between "snapshots", defaults to 1
    """

    def __init__(self, model, cancerCellName="CancerCell", interval=1):

        model.output["cancerCellProperties"]["avgOxygen"] = list()
        model.output["cancerCellProperties"]["minOxygen"] = list()
        model.output["cancerCellProperties"]["maxOxygen"] = list()
        model.output["cancerCellProperties"]["OxygenDistributions"] = list()
        self.cancer_cell_name = cancerCellName
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.oxygen_env_name = model.properties["envNames"]["oxygenEnvName"]
        self.interval = interval

    def step_epilogue(self, model):

        cancerCells = [a for a in model.schedule.agents if (
                a.__class__.__name__ == self.cancer_cell_name and not a.dead)
                       or a.__class__.__name__ == "HealthyCell"]
        if len(cancerCells) > 0:
            coordinates = [a.environment_positions[self.agent_env_name] for a
                           in
                           cancerCells]
            concentrations = [model.environments[self.oxygen_env_name].grid[c]
                              for c in coordinates]
            model.output["cancerCellProperties"]["avgOxygen"].append(
                np.mean(concentrations))
            model.output["cancerCellProperties"]["minOxygen"].append(
                min(concentrations))
            model.output["cancerCellProperties"]["maxOxygen"].append(
                max(concentrations))

            if model.current_epoch % self.interval == 0 or \
                    model.current_epoch \
                    == model.epochs - 1:
                n, bins, patches = plt.hist(concentrations)
                model.output["cancerCellProperties"][
                    "OxygenDistributions"].append(
                    {"n": n, "bins": bins, "epoch": model.current_epoch})
        else:
            model.output["cancerCellProperties"]["avgOxygen"].append(0)
            model.output["cancerCellProperties"]["minOxygen"].append(0)
            model.output["cancerCellProperties"]["maxOxygen"].append(0)
