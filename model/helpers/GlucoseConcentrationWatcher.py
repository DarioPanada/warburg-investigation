'''
Keeps track of average, minimum and maximum glucose concentrations at cancer cells at each epoch
'''
from core.Steppables import Helper
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

class GlucoseConcentrationWatcher(Helper, object):

    '''
   :param model - The model object
   :param cancerCellClassName - Name of the class used to represent cancer cells (defaults to CancerCell)
   :param interval - A value n such that a snapshot of glucose concentrations will be taken every n epochs
   (defaults to 1)
   '''
    def __init__(self, model, cancerCellName="CancerCell", interval=1):

        model.output["cancerCellProperties"]["avgGlucose"] = list()
        model.output["cancerCellProperties"]["minGlucose"] = list()
        model.output["cancerCellProperties"]["maxGlucose"] = list()
        model.output["cancerCellProperties"]["GlucoseDistributions"] = list()
        self.cancerCellName = cancerCellName
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.glucoseEnvName = model.properties["envNames"]["glucoseEnvName"]
        self.interval = interval

    def stepEpilogue(self, model):

        cancerCells = [a for a in model.schedule.agents if a.__class__.__name__ == self.cancerCellName and not a.dead]
        if len(cancerCells) > 0:
            coordinates = [a.environmentPositions[self.agentEnvName] for a in cancerCells]
            concentrations = [model.environments[self.glucoseEnvName].grid[c] for c in coordinates]
            model.output["cancerCellProperties"]["avgGlucose"].append(np.mean(concentrations))
            model.output["cancerCellProperties"]["minGlucose"].append(min(concentrations))
            model.output["cancerCellProperties"]["maxGlucose"].append(max(concentrations))

            if model.currentEpoch % self.interval == 0 or model.currentEpoch == model.epochs - 1:
                n, bins, patches = plt.hist(concentrations)
                model.output["cancerCellProperties"]["GlucoseDistributions"].append(
                    {"n": n, "bins": bins, "epoch": model.currentEpoch})
        else:
            model.output["cancerCellProperties"]["avgGlucose"].append(0)
            model.output["cancerCellProperties"]["minGlucose"].append(0)
            model.output["cancerCellProperties"]["maxGlucose"].append(0)