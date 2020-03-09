import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.switch_backend("agg")

import numpy as np
from panaxea.core.Steppables import Helper


class CancerCellWatcher(Helper, object):
    """
    Cancer cell watcher object. Collects average average hif expression rates,
    vegf secretion rate, metabolic rate
    and pSynthesis. Creates appropriate keys in model.outputs[
    "cancerCellProperties"] to store such values
    at each epoch

    Attributes
    ----------
    model : Model
        The model instance
    cancer_cell_class_name : string, optional
        The name of the cancer cell class, defaults to Cancer Cell
    distribution_interval : int, optional
        A value n such that a snapshot of HIF Rates distribution will be taken
        every n epochs. Optional, defaults to 1
    """

    def __init__(self, model, cancerCellClassName="CancerCell",
                 distributionInterval=1):

        model.output["cancerCellProperties"]["avgHif"] = []
        model.output["cancerCellProperties"]["avgVegf"] = []
        model.output["cancerCellProperties"]["avgMetabolicRates"] = []
        model.output["cancerCellProperties"]["avgPSynthesis"] = []
        model.output["cancerCellProperties"][
            "HIFExpressionRatesDistributions"] = []
        model.output["cancerCellProperties"]["numWarburgCells"] = []

        self.cancer_class_name = cancerCellClassName
        self.distribution_interval = distributionInterval

    def step_epilogue(self, model):

        cancerCells = [a for a in model.schedule.agents if
                       a.__class__.__name__ == self.cancer_class_name and not
                       a.dead]

        if len(cancerCells) == 0:
            print("HERE")
            model.output["cancerCellProperties"]["avgHif"].append(0)
            model.output["cancerCellProperties"]["avgVegf"].append(0)
            model.output["cancerCellProperties"]["avgMetabolicRates"].append(0)
            model.output["cancerCellProperties"]["avgPSynthesis"].append(0)
        else:
            hif_rates = [a.current_hif_rate for a in cancerCells]
            model.output["cancerCellProperties"]["avgHif"].append(
                np.mean(hif_rates))

            vegf_rates = [a.current_vegf_secretion_rate for a in cancerCells]
            model.output["cancerCellProperties"]["avgVegf"].append(
                np.mean(vegf_rates))

            metabolic_rates = [a.current_metabolic_rate for a in cancerCells]
            model.output["cancerCellProperties"]["avgMetabolicRates"].append(
                np.mean(metabolic_rates))

            p_synthesis = [a.current_p_synthesis for a in cancerCells]
            model.output["cancerCellProperties"]["avgPSynthesis"].append(
                np.mean(p_synthesis))

            # Percentage of warburg cells
            warburg_cells = [c for c in cancerCells if c.warburg_switch]
            model.output["cancerCellProperties"]["numWarburgCells"].append(
                float(len(warburg_cells)) / float(len(cancerCells)))

            if model.current_epoch % self.distribution_interval == 0 or \
                    model.current_epoch == model.epochs - 1:
                n, bins, patches = plt.hist(hif_rates,
                                            bins=list(np.arange(0, 17, 1)))
                model.output["cancerCellProperties"][
                    "HIFExpressionRatesDistributions"].append(
                    {"n": n, "bins": bins, "epoch": model.current_epoch})
