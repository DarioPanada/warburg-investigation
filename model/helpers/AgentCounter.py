from panaxea.core.Steppables import Helper


class AgentCounter(Helper, object):
    """
    Counts number of agent types at each epoch and saves them to model output.

    These are saved under the key agentNums as a map of keys, one per class
    name, to lists, where at each epoch the number of class instances is
    appended.

    Attributes
    ----------
    model : Model
        The model instance
    cancer_cell_class_name : string, optional
        The name of the cancer cell class, defaults to Cancer Cell
    """

    def __init__(self, model, cancer_cell_class_name="CancerCell"):
        self.cancer_cell_class_name = cancer_cell_class_name
        self.tip_cell_class_name = "TipCell"
        self.trunk_cell_name = "TrunkCell"

        agent_nums = {
            "cancerCells": [],
            "aliveCancerCells": [],
            "deadCancerCells": [],
            "tipCells": []
        }

        model.output["agentNums"] = agent_nums

    def step_epilogue(self, model):
        cancer_cells = len([a for a in model.schedule.agents if
                            a.__class__.__name__ ==
                            self.cancer_cell_class_name])
        tip_cells = len([a for a in model.schedule.agents if
                         a.__class__.__name__ == self.tip_cell_class_name or
                         a.__class__.__name__ == self.trunk_cell_name])
        alive_cancer_cells = len([a for a in model.schedule.agents if
                                  a.__class__.__name__ ==
                                  self.cancer_cell_class_name and not a.dead])
        dead_cancer_cells = len([a for a in model.schedule.agents if
                                 a.__class__.__name__ ==
                                 self.cancer_cell_class_name and a.dead])

        model.output["agentNums"]["cancerCells"].append(cancer_cells)
        model.output["agentNums"]["tipCells"].append(tip_cells)
        model.output["agentNums"]["deadCancerCells"].append(dead_cancer_cells)
        model.output["agentNums"]["aliveCancerCells"].append(
            alive_cancer_cells)
