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
        self.update_oxygen_emission_rate()

        self.base_glucose_secretion_rate = \
            model.properties["agents"]["endothelialCells"][
                "glucoseSecretionRate"]
        self.update_glucose_secretion_rate()

    def step_main(self, model):
        self.update_oxygen_emission_rate()

    def update_glucose_secretion_rate(self):
        self.glucose_secretion_rate = self.radius * \
                                      self.base_glucose_secretion_rate

    def update_oxygen_emission_rate(self):
        self.oxygen_emission_rate = self.radius * self.baseOxygenEmissionRate


class TipCell(EndothelialCell, object):
    def __init__(self, model, radius=1):
        super(TipCell, self).__init__(model, radius=radius)

        # Copying environment names to avoid continuosly accessing the model
        self.vegf_env = model.properties["envNames"]["vegfEnvName"]
        self.agent_env = model.properties["envNames"]["agentEnvName"]

        endothelial_cells_properties = model.properties["agents"][
            "endothelialCells"]
        self.minimum_vegf_concentration = endothelial_cells_properties[
            "minimumVegfConcentration"]
        self.division_delay = endothelial_cells_properties["divisionDelay"]

        self.cell_age = 0

    def step_main(self, model):
        super(TipCell, self).step_main(model)
        self.cell_age = min(self.division_delay, self.cell_age + 1)

        current_position = self.environment_positions[self.agent_env]

        taf_at_current_position = model.environments[self.vegf_env].grid[
            current_position]

        if taf_at_current_position >= self.minimum_vegf_concentration and \
                self.cell_age == self.division_delay and 10. * random() < \
                taf_at_current_position:
            self.cell_age = 0

            neigh = model.environments[self.agent_env].get_moore_neighbourhood(
                self.environment_positions[self.agent_env])

            # get taf concentration at each neighbouring position
            neigh_to_taf = [(n, model.environments[self.vegf_env].grid[n]) for
                            n
                            in neigh]

            def f(x):
                return x[1]

            sorted_neigh = sorted(neigh_to_taf, key=f, reverse=True)

            new_pos = sorted_neigh.pop(0)[0]

            while model.environments[self.agent_env].grid[new_pos].__len__() \
                    >= \
                    model.properties["maxAgentDensity"]:
                if len(sorted_neigh) == 0:
                    return
                else:
                    new_pos = sorted_neigh.pop(0)[0]

            self.move_agent(self.agent_env, new_pos, model)

            # Create tip cell at old position1
            t = TrunkCell(model, radius=self.radius)
            t.add_agent_to_grid(self.agent_env, current_position, model)
            model.schedule.agents_to_schedule.add(t)

    # A cell automatically sprouts if its radius is > 1
    def _decide_sprout_linear(self, tafConcentration):

        return self.radius > 1

    # Given a list of sorted neighbours where each takes the form (coords,
    # score), returns a neighbour
    # where those with a higher score have the highest probability of being
    # picked
    def _get_next_neigh_from_scored(self, neighs):
        threshold = random()

        for n in neighs:
            if n[1] > threshold:
                return n[0]

    # Assigns a score to each neighbourhood position where 0 = awful and 1 =
    # best and sum(scores) = 1
    def _rank_neighbours(self, model):
        neigh = model.environments[self.agent_env].get_moore_neighbourhood(
            self.environment_positions[self.agent_env])

        # get taf concentration at each neighbouring position
        neigh_to_taf = [(n, model.environments[self.vegf_env].grid[n]) for n in
                        neigh]

        # get total taf concentration
        total_taf = sum([x[1] for x in neigh_to_taf])

        # In the case where there is no taf, we treat all neighbouring
        # positions alike
        if total_taf < 0.0001:
            scored_neigh = [(x[0], 0) for x in neigh_to_taf]
        else:
            # scoring the neighbourhood
            scored_neigh = [(x[0], float(x[1]) / float(total_taf)) for x in
                            neigh_to_taf]

        # sorting by taf concentration
        def f(x):
            return x[1]

        sorted_neigh = sorted(scored_neigh, key=f)

        accumulator = 0
        accumulated_neigh = []

        for n in sorted_neigh:
            accumulator = accumulator + n[1]
            accumulated_neigh.append((n[0], accumulator))

        return accumulated_neigh


class TrunkCell(EndothelialCell, object):
    def __init__(self, model, radius=1):
        super(TrunkCell, self).__init__(model, radius=radius)
