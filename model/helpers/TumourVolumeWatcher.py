from cmath import sqrt

from panaxea.core.Steppables import Helper


class TumourVolumeWatcher(Helper, object):
    """
    At each time-step, calculates an estimate of the tumour volume by
    looking at
    the euclidean distance of the two farthest
    cancer cells.
    """

    def __init__(self, model, cancer_cell_class_name="CancerCell"):
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.cancer_cell_class_name = cancer_cell_class_name

        model.output["maxDistances"] = []

    def step_epilogue(self, model):
        cancer_cells_coords = [a.environment_positions[self.agent_env_name]
                               for a in model.schedule.agents
                               if a.__class__.__name__ ==
                               self.cancer_cell_class_name]
        scored_coords = [(sum(c), c) for c in cancer_cells_coords]

        def f(c):
            return c[0]

        if len(scored_coords) > 0:
            max_coord = max(scored_coords, key=f)[1]
            min_coord = min(scored_coords, key=f)[1]

            dist = sqrt((max_coord[0] - min_coord[0]) ** 2 + (
                    max_coord[1] - min_coord[1]) ** 2 + (
                                max_coord[2] - min_coord[2]) ** 2)
        else:
            dist = 0

        model.output["maxDistances"].append(dist)
