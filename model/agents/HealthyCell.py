from core.Steppables import Agent


class HealthyCell(Agent):

    def __init__(self, model):
        super(HealthyCell, self).__init__()
        self.currentMetabolicRate = model.properties["agents"]["healthyTissues"]["oxygenUptakeRate"]
        self.glucoseUptakeRate =  model.properties["agents"]["cancerCells"]["minGlucoseUptakeRate"]
        self.dead = False
