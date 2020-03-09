import random
from numpy.polynomial import Polynomial
from panaxea.core.Steppables import Agent


class CancerCell(Agent, object):
    def __init__(self, model, warburgSwitch=False):
        super(CancerCell, self).__init__()

        self.current_state = "G1"
        self.progress_in_state = 0
        self.cell_cycle_length = model.properties["agents"][
            "baseCellCycleLength"]
        self.cell_cycle_order = ["G1", "S", "G2", "M"]
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.oxygen_env_name = model.properties["envNames"]["oxygenEnvName"]
        self.glucose_env_name = model.properties["envNames"]["glucoseEnvName"]
        self.drug_env_name = model.properties["envNames"]["drugEnvName"]

        cancer_cell_props = model.properties["agents"]["cancerCells"]

        self.base_hif_rate = cancer_cell_props["baseHifRate"]
        self.oxygen_to_hif_coeffs_hypoxic = \
            cancer_cell_props["oxygenToHifCoeffs"]["hypoxic"]
        self.oxygen_to_hif_coeffs_ultra_hypoxic = \
            cancer_cell_props["oxygenToHifCoeffs"]["ultraHypoxic"]
        self.oxygen_to_hif_coeffs_hypoxic_warburg = \
            cancer_cell_props["oxygenToHifCoeffs"]["warburg"]
        self.oxygen_hypoxic_domain = cancer_cell_props["domains"]["hypoxic"]
        self.oxygen_ultra_hypoxic_domain = cancer_cell_props["domains"][
            "ultraHypoxic"]
        self.oxygen_warburg_hypoxic_domain = cancer_cell_props["domains"][
            "warburgHypoxic"]

        self.current_hif_rate = cancer_cell_props["minHIF"]

        self.hif_to_metabolic_rate_coeffs = cancer_cell_props[
            "hifToMetabolicRateCoeffs"]
        self.current_metabolic_rate = 1

        self.hif_to_proliferation_rate_coeffs = cancer_cell_props[
            "hifToProliferationRateCoeffs"]
        self.min_p_synthesis = cancer_cell_props["minPSynthesis"]
        self.current_p_synthesis = self.min_p_synthesis

        self.hif_to_vegf_secretion_rate_coeffs = cancer_cell_props[
            "hifToVegfSecretionRateCoeffs"]
        self.current_vegf_secretion_rate = 1

        self.minimum_oxygen_concentration = cancer_cell_props[
            "minimumOxygenConcentration"]
        self.min_glucose_warburg = cancer_cell_props["minGlucoseWarburg"]
        self.min_glucose_non_warburg = cancer_cell_props[
            "minGlucoseNonWarburg"]

        self.hif_range = cancer_cell_props["HIFRange"]

        self.min_glucose_uptake_rate = cancer_cell_props[
            "minGlucoseUptakeRate"]
        self.max_glucose_uptake_rate = cancer_cell_props[
            "maxGlucoseUptakeRate"]
        self.glucose_uptake_rate = self.min_glucose_uptake_rate

        self.p_warburg_switch = cancer_cell_props["pWarburgSwitch"]
        self.warburg_switch = warburgSwitch

        self.dead = False
        self.quiescent = False

        self.drug_at_pos = 0

        # To allow for different drug effects to be tested, this will be
        # loaded as part of the config.
        self.max_vegf = cancer_cell_props["maxVegfSecretionRate"]
        self.age = 0
        self.min_hif = cancer_cell_props["minHIF"]

    def progress_cell_(self, model):
        # time to divide
        if self.current_state == self.cell_cycle_order[-1] and \
                self.progress_in_state == \
                self.cell_cycle_length[self.current_state]:
            current_pos = self.environment_positions[self.agent_env_name]
            target_pos = None

            if len(model.environments[self.agent_env_name].grid[
                       (current_pos[0], current_pos[1], current_pos[2])]) < \
                    model.properties["maxAgentDensity"]:
                target_pos = current_pos
            else:
                moore_target = model.environments[
                    self.agent_env_name].getLeastPopulatedMooreNeigh(
                    current_pos)
                if len(model.environments[self.agent_env_name].grid[
                           (moore_target[0], moore_target[1], moore_target[
                               2])]) < \
                        model.properties["maxAgentDensity"]:
                    target_pos = moore_target

            if target_pos is not None:
                # Creating new cancer cell and adding it at current position
                # in the designated environment
                c = CancerCell(model)
                c.add_agent_to_grid(self.agent_env_name, target_pos, model)
                model.schedule.agents_to_schedule.add(c)

                # Resetting current cell
                self.progress_in_state = 0
                self.current_state = self.cell_cycle_order[0]
                self.warburg_switch = False
                self.age = 0

        # progress in cell life-cyle
        elif self.progress_in_state == \
                self.cell_cycle_length[self.current_state]:
            if self.current_state == "G1" and random.random() > \
                    self.current_p_synthesis:
                pass
            else:
                self.progress_in_state = 0
                # Get index of current state and set next state
                current_index = self.cell_cycle_order.index(self.current_state)
                self.current_state = self.cell_cycle_order[current_index + 1]
        else:
            self.progress_in_state = self.progress_in_state + 1

    def step_main(self, model):
        current_pos = self.environment_positions[self.agent_env_name]
        self.oxygen_at_pos = model.environments[self.oxygen_env_name].grid[
            current_pos]
        self.glucose_at_pos = model.environments[self.glucose_env_name].grid[
            current_pos]
        self.drug_at_pos = model.environments[self.drug_env_name].grid[
            current_pos]

        # No logic executed on dead or quiescent cells
        if not self.dead and not self.quiescent:
            self.age = self.age + 1

            if random.random() < self.p_warburg_switch and not \
                    self.warburg_switch:
                self.warburg_switch = True
                self.glucose_uptake_rate = self.max_glucose_uptake_rate
            # If the cell is not warburg, then in order for it to survive it
            # must:
            # Have access to minimum oxygen concentration
            # Have access to minimum glucose for non warburg cells
            # If the cell is a warburg cell then, in order for it to
            # survive, it must:
            # Have access to minimum glucose for warburg cells
            if not self.decide_die_():

                # If a cell is a warburg cell and the dosage at position is
                # greater than the minimum threshold,
                # it becomes quiescent. (We return so that no further cell
                # logic is executed)

                # TEMPORARILY DISABLING DRUGS, THIS SHOULD BE TURNED INTO A
                # FLAG THAT CAN BE SET FROM PARAMS

                self._update_hif_and_mediated(model)
                self.progress_cell_(model)
            else:
                # A dead cell is not quiescent
                self.dead = True
                self.quiescent = False

    def decide_die_(self):
        if self.warburg_switch and self.glucose_at_pos < \
                self.min_glucose_warburg:
            self.cause_of_death = {
                "cause": "Lack of glucose",
                "glucoseAtPos": self.glucose_at_pos,
                "warburg": self.warburg_switch,
                "age": self.age
            }
            return True
        elif not self.warburg_switch:
            if self.oxygen_at_pos < self.minimum_oxygen_concentration:
                self.cause_of_death = {
                    "cause": "Lack of oxygen",
                    "oxygenAtPos": self.oxygen_at_pos,
                    "warburg": self.warburg_switch,
                    "age": self.age
                }
                return True
            elif self.glucose_at_pos < self.min_glucose_non_warburg:
                self.cause_of_death = {
                    "cause": "Lack of glucose",
                    "glucoseAtPos": self.glucose_at_pos,
                    "warburg": self.warburg_switch,
                    "age": self.age
                }
                return True

        return False

    def _update_hif_and_mediated(self, model):
        self._update_hif_expression_rate(model)
        self._update_metabolic_rate()
        self._update_p_synthesis()
        self._update_vegf_secretion_rate()

    def _update_p_synthesis(self):

        p = Polynomial(coef=self.hif_to_proliferation_rate_coeffs,
                       domain=self.hif_range)

        p_synthesis = p(self.current_hif_rate)

        self.current_p_synthesis = p_synthesis

    def _update_vegf_secretion_rate(self):

        p = Polynomial(coef=self.hif_to_vegf_secretion_rate_coeffs,
                       domain=self.hif_range)

        vegf_rate = p(self.current_hif_rate)

        adjusted_vegf_rate = max(0, min(1, vegf_rate))
        self.current_vegf_secretion_rate = self.max_vegf * adjusted_vegf_rate

    def _update_metabolic_rate(self):

        p = Polynomial(coef=self.hif_to_metabolic_rate_coeffs,
                       domain=self.hif_range)

        metabolic_rate = p(self.current_hif_rate)

        self.current_metabolic_rate = metabolic_rate

    def _update_hif_expression_rate(self, model):

        if not self.warburg_switch:
            new_rate = self._calculate_hif_expression_rate_from_oxygen(
                self.oxygen_at_pos)
        else:
            new_rate = \
                self._calculate_hif_expression_rate_from_oxygen_warburg(
                    self.oxygen_at_pos)

        max_shift_per_epoch = 0.1

        current_rate = self.current_hif_rate

        if current_rate > new_rate:
            new_rate = max(0, current_rate - min(max_shift_per_epoch,
                                                 current_rate - new_rate))
        else:
            new_rate = min(
                current_rate + min(max_shift_per_epoch,
                                   new_rate - current_rate), 16)

        self.current_hif_rate = self.base_hif_rate * new_rate

    def _calculate_hif_expression_rate_from_oxygen(self, oxygen_at_pos):

        if oxygen_at_pos > self.oxygen_hypoxic_domain:
            return 1

        if oxygen_at_pos > self.oxygen_ultra_hypoxic_domain:
            coef = self.oxygen_to_hif_coeffs_hypoxic
            domain = [self.oxygen_ultra_hypoxic_domain,
                      self.oxygen_hypoxic_domain]
        else:
            coef = self.oxygen_to_hif_coeffs_ultra_hypoxic
            domain = [0.0, self.oxygen_ultra_hypoxic_domain]

        p = Polynomial(coef=coef, domain=domain)
        return p(oxygen_at_pos)

    def _calculate_hif_expression_rate_from_oxygen_warburg(self,
                                                           oxygen_at_pos):
        if oxygen_at_pos > self.oxygen_warburg_hypoxic_domain:
            return self.min_hif

        if oxygen_at_pos > self.oxygen_ultra_hypoxic_domain:
            coef = self.oxygen_to_hif_coeffs_hypoxic_warburg
            domain = [self.oxygen_ultra_hypoxic_domain,
                      self.oxygen_warburg_hypoxic_domain]
        else:
            coef = self.oxygen_to_hif_coeffs_ultra_hypoxic
            domain = [0.0, self.oxygen_ultra_hypoxic_domain]

        p = Polynomial(coef=coef, domain=domain)
        return p(oxygen_at_pos)
