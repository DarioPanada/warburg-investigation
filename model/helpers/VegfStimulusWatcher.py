'''
Keeps track of average vegf concentration at blood vessels
'''
import numpy as np
from core.Steppables import Helper


class VegfStimulusWatcher(Helper, object):

    def __init__(self, model, endothelialCellNames=["TipCell"]):
        model.output["endothelialCellProperties"] = dict()
        model.output["endothelialCellProperties"]["avgVegf"] = []
        self.endothelialCellNames = endothelialCellNames
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.vegfEnvName = model.properties["envNames"]["vegfEnvName"]

    def step_epilogue(self, model):
        edothelialCells = [a for a in model.schedule.agents if
                           a.__class__.__name__ in self.endothelialCellNames]
        coordinates = [a.environment_positions[self.agentEnvName] for a in
                       edothelialCells]
        concentrations = [model.environments[self.vegfEnvName].grid[c] for c in
                          coordinates]
        model.output["endothelialCellProperties"]["avgVegf"].append(
            np.mean(concentrations))
