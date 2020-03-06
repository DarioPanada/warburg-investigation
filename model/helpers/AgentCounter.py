from core.Steppables import Helper

'''
Counts number of agent types at each epoch and saves them to model output
'''
class AgentCounter(Helper, object):

    def __init__(self, model, cancerCellClassName="CancerCell"):
        self.cancerCellClassName = cancerCellClassName
        self.tipCellClassName = "TipCell"
        self.trunkCellName = "TrunkCell"

        agentNums = {
            "cancerCells": [],
            "aliveCancerCells" : [],
            "deadCancerCells" : [],
            "tipCells": []
        }

        model.output["agentNums"] = agentNums

    def stepEpilogue(self, model):
        cancerCells = len([a for a in model.schedule.agents if a.__class__.__name__ == self.cancerCellClassName])
        tipCells = len([a for a in model.schedule.agents if a.__class__.__name__ == self.tipCellClassName or a.__class__.__name__ == self.trunkCellName])
        aliveCancerCells = len([a for a in model.schedule.agents if a.__class__.__name__ == self.cancerCellClassName and not a.dead])
        deadCancerCells = len([a for a in model.schedule.agents if a.__class__.__name__ == self.cancerCellClassName and a.dead])


        model.output["agentNums"]["cancerCells"].append(cancerCells)
        model.output["agentNums"]["tipCells"].append(tipCells)
        model.output["agentNums"]["deadCancerCells"].append(deadCancerCells)
        model.output["agentNums"]["aliveCancerCells"].append(aliveCancerCells)