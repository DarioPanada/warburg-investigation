'''
Keeps track of average, minimum and maximum oxygen concentration at cancer
cells at each epoch
'''
import matplotlib
import numpy as np
from panaxea.core.Steppables import Helper

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class OxygenConcentrationWatcher(Helper, object):
    '''
    :param model - The model object
    :param cancerCellClassName - Name of the class used to represent cancer
    cells (defaults to CancerCell)
    :param interval - A value n such that a snapshot of oxygen
    concentrations will be taken every n epochs
    (defaults to 1)
    '''

    def __init__(self, model, cancerCellName="CancerCell", interval=1):

        model.output["cancerCellProperties"]["avgOxygen"] = list()
        model.output["cancerCellProperties"]["minOxygen"] = list()
        model.output["cancerCellProperties"]["maxOxygen"] = list()
        model.output["cancerCellProperties"]["OxygenDistributions"] = list()
        self.cancerCellName = cancerCellName
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.oxygenEnvName = model.properties["envNames"]["oxygenEnvName"]
        self.interval = interval

    def step_epilogue(self, model):

        cancerCells = [a for a in model.schedule.agents if (
                    a.__class__.__name__ == self.cancerCellName and not a.dead)
                       or a.__class__.__name__ == "HealthyCell"]
        if len(cancerCells) > 0:
            coordinates = [a.environment_positions[self.agentEnvName] for a in
                           cancerCells]
            concentrations = [model.environments[self.oxygenEnvName].grid[c]
                              for c in coordinates]
            model.output["cancerCellProperties"]["avgOxygen"].append(
                np.mean(concentrations))
            model.output["cancerCellProperties"]["minOxygen"].append(
                min(concentrations))
            model.output["cancerCellProperties"]["maxOxygen"].append(
                max(concentrations))

            if model.current_epoch % self.interval == 0 or model.current_epoch\
                    == model.epochs - 1:
                n, bins, patches = plt.hist(concentrations)
                model.output["cancerCellProperties"][
                    "OxygenDistributions"].append(
                    {"n": n, "bins": bins, "epoch": model.current_epoch})
        else:
            model.output["cancerCellProperties"]["avgOxygen"].append(0)
            model.output["cancerCellProperties"]["minOxygen"].append(0)
            model.output["cancerCellProperties"]["maxOxygen"].append(0)
