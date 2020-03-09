from panaxea.core.Steppables import Agent


class HealthyCell(Agent):

    def __init__(self, model):
        super(HealthyCell, self).__init__()
        self.current_metabolic_rate = \
            model.properties["agents"]["healthyTissues"]["oxygenUptakeRate"]
        self.glucose_uptake_rate = model.properties["agents"]["cancerCells"][
            "minGlucoseUptakeRate"]
        self.dead = False
