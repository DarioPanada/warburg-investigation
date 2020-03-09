import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.switch_backend("agg")
import numpy as np
from panaxea.core.Steppables import Helper


class GlucoseConcentrationWatcher(Helper, object):
    """
    Keeps track of average, minimum and maximum glucose concentrations at
    cancer cells at each epoch.

    Values are stored in the output under the key cancerCellProperties.

    Attributes
    ----------
    model : Model
        The current model instance
    cancer_cell_class_name : string, optional
        The name of the cancer cell class, defaults to Cancer Cell
        distribution_interval : int, optional
    interval : int, optional
        A value n such that a snapshot of glucose
        concentrations will be taken every n epochs. Optional, defaults to 1
    """

    def __init__(self, model, cancer_cell_name="CancerCell", interval=1):

        model.output["cancerCellProperties"]["avgGlucose"] = list()
        model.output["cancerCellProperties"]["minGlucose"] = list()
        model.output["cancerCellProperties"]["maxGlucose"] = list()
        model.output["cancerCellProperties"]["GlucoseDistributions"] = list()
        self.cancer_cell_name = cancer_cell_name
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.glucose_env_name = model.properties["envNames"]["glucoseEnvName"]
        self.interval = interval

    def step_epilogue(self, model):

        cancer_cells = [a for a in model.schedule.agents if
                        a.__class__.__name__ == self.cancer_cell_name and not
                        a.dead]
        if len(cancer_cells) > 0:
            coordinates = [a.environment_positions[self.agent_env_name] for a
                           in
                           cancer_cells]
            concentrations = [model.environments[self.glucose_env_name].grid[c]
                              for c in coordinates]
            model.output["cancerCellProperties"]["avgGlucose"].append(
                np.mean(concentrations))
            model.output["cancerCellProperties"]["minGlucose"].append(
                min(concentrations))
            model.output["cancerCellProperties"]["maxGlucose"].append(
                max(concentrations))

            if model.current_epoch % self.interval == 0 or \
                    model.current_epoch \
                    == model.epochs - 1:
                n, bins, patches = plt.hist(concentrations)
                model.output["cancerCellProperties"][
                    "GlucoseDistributions"].append(
                    {"n": n, "bins": bins, "epoch": model.current_epoch})
        else:
            model.output["cancerCellProperties"]["avgGlucose"].append(0)
            model.output["cancerCellProperties"]["minGlucose"].append(0)
            model.output["cancerCellProperties"]["maxGlucose"].append(0)
