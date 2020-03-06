'''
Keeps track of average vegf concentration at blood vessels
'''
from core.Steppables import Helper
import numpy as np

class VegfStimulusWatcher(Helper, object):

    def __init__(self, model, endothelialCellNames = ["TipCell"]):

        model.output["endothelialCellProperties"] = dict()
        model.output["endothelialCellProperties"]["avgVegf"] = []
        self.endothelialCellNames = endothelialCellNames
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.vegfEnvName = model.properties["envNames"]["vegfEnvName"]

    def stepEpilogue(self, model):

        edothelialCells = [a for a in model.schedule.agents if a.__class__.__name__ in self.endothelialCellNames]
        coordinates = [a.environmentPositions[self.agentEnvName] for a in edothelialCells]
        concentrations = [model.environments[self.vegfEnvName].grid[c] for c in coordinates]
        model.output["endothelialCellProperties"]["avgVegf"].append(np.mean(concentrations))