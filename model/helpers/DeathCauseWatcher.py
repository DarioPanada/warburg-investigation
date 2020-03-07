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
            cancerCells = [a for a in model.schedule.agents if
                           a.__class__.__name__ == "CancerCell" and a.dead]

            warburgDeathGlucose = [c for c in cancerCells if
                                   c.warburgSwitch and c.causeOfDeath[
                                       "cause"] == "Lack of glucose"]
            warburgDeathOxygen = [c for c in cancerCells if
                                  c.warburgSwitch and c.causeOfDeath[
                                      "cause"] == "Lack of oxygen"]
            nonWarburgDeathGlucose = [c for c in cancerCells if
                                      not c.warburgSwitch and c.causeOfDeath[
                                          "cause"] == "Lack of glucose"]
            nonWarburgDeathOxygen = [c for c in cancerCells if
                                     not c.warburgSwitch and c.causeOfDeath[
                                         "cause"] == "Lack of oxygen"]

            if (len(warburgDeathGlucose) + len(warburgDeathOxygen) + len(
                    nonWarburgDeathGlucose) + len(
                    nonWarburgDeathOxygen)) != len(cancerCells):
                print(
                    "ERROR, LENGTH OF COLLECTED DEAD CELLS DIFFERENT FROM "
                    "LENGTH OF TOTAL DEAD CELLS")
                exit()

            # summarising
            summary = {
                "warburgDeathGlucose": {"num": len(warburgDeathGlucose),
                                        "avgAge": np.mean([c.age for c in
                                                           warburgDeathGlucose]),
                                        "stDev": np.std([c.age for c in
                                                         warburgDeathGlucose])},
                "warburgDeathOxygen": {"num": len(warburgDeathOxygen),
                                       "avgAge": np.mean([c.age for c in
                                                          warburgDeathOxygen]),
                                       "stDev": np.std([c.age for c in
                                                        warburgDeathOxygen])},
                "nonWarburgDeathOxygen": {"num": len(nonWarburgDeathOxygen),
                                          "avgAge": np.mean([c.age for c in
                                                             nonWarburgDeathOxygen]),
                                          "stDev": np.std([c.age for c in
                                                           nonWarburgDeathOxygen])},
                "nonWarburgDeathGlucose": {"num": len(nonWarburgDeathGlucose),
                                           "avgAge": np.mean([c.age for c in
                                                              nonWarburgDeathGlucose]),
                                           "stDev": np.std([c.age for c in
                                                            nonWarburgDeathGlucose])},

            }

            model.output["causesOfDeath"].append(summary)
