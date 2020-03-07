from panaxea.core.Steppables import Agent
from random import random


class EndothelialCell(Agent, object):

    def __init__(self, model, radius=1):
        super(EndothelialCell, self).__init__()

        if radius < 1:
            raise Exception(
                'Trying to instanstiate endothelial cell with invalid radius '
                '%s' % str(
                    radius))

        self.radius = radius
        self.baseOxygenEmissionRate = \
        model.properties["agents"]["endothelialCells"][
            "baseOxygenEmissionRate"]
        self.updateOxygenEmissionRate()

        self.baseGlucoseSecretionRate = \
        model.properties["agents"]["endothelialCells"]["glucoseSecretionRate"]
        self.updateGlucoseSecretionRate()

    def stepMain(self, model):
        self.updateOxygenEmissionRate()

    def updateGlucoseSecretionRate(self):
        self.glucoseSecretionRate = self.radius * self.baseGlucoseSecretionRate

    def updateOxygenEmissionRate(self):
        self.oxygenEmissionRate = self.radius * self.baseOxygenEmissionRate


class TipCell(EndothelialCell, object):
    def __init__(self, model, radius=1):
        super(TipCell, self).__init__(model, radius=radius)

        # Copying environment names to avoid continuosly accessing the model
        self.vegfEnv = model.properties["envNames"]["vegfEnvName"]
        self.agentEnv = model.properties["envNames"]["agentEnvName"]

        endothelialCellsProperties = model.properties["agents"][
            "endothelialCells"]
        self.minimumVegfConcentration = endothelialCellsProperties[
            "minimumVegfConcentration"]
        self.divisionDelay = endothelialCellsProperties["divisionDelay"]

        self.cellAge = 0

    def step_main(self, model):
        super(TipCell, self).stepMain(model)
        self.cellAge = min(self.divisionDelay, self.cellAge + 1)

        currentPosition = self.environment_positions[self.agentEnv]

        tafAtCurrentPosition = model.environments[self.vegfEnv].grid[
            currentPosition]

        if tafAtCurrentPosition >= self.minimumVegfConcentration and \
                self.cellAge == self.divisionDelay and 10. * random() < \
                tafAtCurrentPosition:
            self.cellAge = 0

            neigh = model.environments[self.agentEnv].getMooreNeighbourhood(
                self.environment_positions[self.agentEnv])

            # get taf concentration at each neighbouring position
            neighToTaf = [(n, model.environments[self.vegfEnv].grid[n]) for n
                          in neigh]

            def f(x):
                return x[1]

            sortedNeigh = sorted(neighToTaf, key=f, reverse=True)

            newPos = sortedNeigh.pop(0)[0]

            while model.environments[self.agentEnv].grid[newPos].__len__() >= \
                    model.properties["maxAgentDensity"]:
                if len(sortedNeigh) == 0:
                    return
                else:
                    newPos = sortedNeigh.pop(0)[0]

            self.MoveAgent(self.agentEnv, newPos, model)

            # Create tip cell at old position1
            t = TrunkCell(model, radius=self.radius)
            t.add_agent_to_grid(self.agentEnv, currentPosition, model)
            model.schedule.agentsToSchedule.add(t)

    # A cell automatically sprouts if its radius is > 1
    def decideSproutLinear_(self, tafConcentration):

        return self.radius > 1

    # Given a list of sorted neighbours where each takes the form (coords,
    # score), returns a neighbour
    # where those with a higher score have the highest probability of being
    # picked
    def getNextNeighFromScored_(self, neighs):
        threshold = random()

        for n in neighs:
            if n[1] > threshold:
                return n[0]

    # Assigns a score to each neighbourhood position where 0 = awful and 1 =
    # best and sum(scores) = 1
    def rankNeighbours_(self, model):
        neigh = model.environments[self.agentEnv].getMooreNeighbourhood(
            self.environment_positions[self.agentEnv])

        # get taf concentration at each neighbouring position
        neighToTaf = [(n, model.environments[self.vegfEnv].grid[n]) for n in
                      neigh]

        # get total taf concentration
        totalTaf = sum([x[1] for x in neighToTaf])

        # In the case where there is no taf, we treat all neighbouring
        # positions alike
        if totalTaf < 0.0001:
            scoredNeigh = [(x[0], 0) for x in neighToTaf]
        else:
            # scoring the neighbourhood
            scoredNeigh = [(x[0], float(x[1]) / float(totalTaf)) for x in
                           neighToTaf]

        # sorting by taf concentration
        def f(x):
            return x[1]

        sortedNeigh = sorted(scoredNeigh, key=f)

        accumulator = 0
        accumulatedNeigh = []

        for n in sortedNeigh:
            accumulator = accumulator + n[1]
            accumulatedNeigh.append((n[0], accumulator))

        return accumulatedNeigh


class TrunkCell(EndothelialCell, object):
    def __init__(self, model, radius=1):
        super(TrunkCell, self).__init__(model, radius=radius)
