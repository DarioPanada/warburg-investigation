import numpy as np
from panaxea.core.Steppables import Helper


class VegfStimulusWatcher(Helper, object):
    '''
    Keeps track of average vegf concentration at blood vessels
    '''

    def __init__(self, model, endothelialCellNames=["TipCell"]):
        model.output["endothelialCellProperties"] = dict()
        model.output["endothelialCellProperties"]["avgVegf"] = []
        self.endothelial_cell_names = endothelialCellNames
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.vegf_env_name = model.properties["envNames"]["vegfEnvName"]

    def step_epilogue(self, model):
        edothelial_cells = [a for a in model.schedule.agents
                            if a.__class__.__name__ in
                            self.endothelial_cell_names]
        coordinates = [a.environment_positions[self.agent_env_name] for a in
                       edothelial_cells]
        concentrations = [model.environments[self.vegf_env_name].grid[c]
                          for c in coordinates]
        model.output["endothelialCellProperties"]["avgVegf"].append(
            np.mean(concentrations))
