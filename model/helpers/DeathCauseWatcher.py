import numpy as np
from panaxea.core.Steppables import Helper

'''
Used to track at each epoch what the CUMULATIVE cause of cancer cell death is.
'''


class DeathCauseWatcher(Helper, object):

    def __init__(self, model, interval):
        model.output["causesOfDeath"] = []
        self.interval = interval

    def step_epilogue(self, model):

        if model.current_epoch % self.interval == 0:
            cancer_cells = [a for a in model.schedule.agents if
                            a.__class__.__name__ == "CancerCell" and a.dead]

            warburg_death_glucose = [c for c in cancer_cells if
                                     c.warburg_switchand and c.cause_of_death[
                                         "cause"] == "Lack of glucose"]
            warburg_death_oxygen = [c for c in cancer_cells if
                                    c.warburg_switch and c.cause_of_death[
                                        "cause"] == "Lack of oxygen"]
            non_warburg_death_glucose = [c for c in cancer_cells if
                                         not c.warburg_switch and
                                         c.cause_of_death[
                                             "cause"] == "Lack of glucose"]
            non_warburg_death_oxygen = [c for c in cancer_cells if
                                        not c.warburg_switch and
                                        c.cause_of_death[
                                            "cause"] == "Lack of oxygen"]

            sum_partial_lengths = sum(
                [len(warburg_death_glucose), len(warburg_death_oxygen),
                 len(non_warburg_death_glucose),
                 len(non_warburg_death_oxygen)]
            )

            if sum_partial_lengths != len(cancer_cells):
                print(
                    "ERROR, LENGTH OF COLLECTED DEAD CELLS DIFFERENT FROM "
                    "LENGTH OF TOTAL DEAD CELLS")
                exit()

            # summarising
            summary = {
                "warburgDeathGlucose": {
                    "num": len(warburg_death_glucose),
                    "avgAge": np.mean([c.age for c in warburg_death_glucose]),
                    "stDev": np.std([c.age for c in warburg_death_glucose])},
                "warburgDeathOxygen": {
                    "num": len(warburg_death_oxygen),
                    "avgAge": np.mean([c.age for c in warburg_death_oxygen]),
                    "stDev": np.std([c.age for c in warburg_death_oxygen])},
                "nonWarburgDeathOxygen": {
                    "num": len(non_warburg_death_oxygen),
                    "avgAge": np.mean(
                        [c.age for c in non_warburg_death_oxygen]),
                    "stDev": np.std(
                        [c.age for c in non_warburg_death_oxygen])},
                "nonWarburgDeathGlucose": {
                    "num": len(non_warburg_death_glucose),
                    "avgAge": np.mean(
                        [c.age for c in non_warburg_death_glucose]),
                    "stDev": np.std([c.age for c in non_warburg_death_glucose])
                },

            }

            model.output["causesOfDeath"].append(summary)
