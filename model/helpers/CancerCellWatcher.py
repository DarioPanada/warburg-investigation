import matplotlib
import numpy as np
from panaxea.core.Steppables import Helper

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.switch_backend("agg")
'''
Cancer cell watcher object. Collects average average hif expression rates, 
vegf secretion rate, metabolic rate
and pSynthesis. Creates appropriate keys in model.outputs[
"cancerCellProperties"] to store such values
at each epoch
'''


class CancerCellWatcher(Helper, object):
    '''
    :param model - The model object
    :param cancerCellClassName - Name of the class used to represent cancer
    cells (defaults to CancerCell)
    :param distributionInterval - A value n such that a snapshot of HIF
    Rates distribution will be taken every n epochs
    (defaults to 1)
    '''

    def __init__(self, model, cancerCellClassName="CancerCell",
                 distributionInterval=1):

        model.output["cancerCellProperties"]["avgHif"] = []
        model.output["cancerCellProperties"]["avgVegf"] = []
        model.output["cancerCellProperties"]["avgMetabolicRates"] = []
        model.output["cancerCellProperties"]["avgPSynthesis"] = []
        model.output["cancerCellProperties"][
            "HIFExpressionRatesDistributions"] = []
        model.output["cancerCellProperties"]["numWarburgCells"] = []

        self.cancerClassName = cancerCellClassName
        self.distributionInterval = distributionInterval

    def step_epilogue(self, model):

        cancerCells = [a for a in model.schedule.agents if
                       a.__class__.__name__ == self.cancerClassName and not
                       a.dead]

        if len(cancerCells) == 0:
            print("HERE")
            model.output["cancerCellProperties"]["avgHif"].append(0)
            model.output["cancerCellProperties"]["avgVegf"].append(0)
            model.output["cancerCellProperties"]["avgMetabolicRates"].append(0)
            model.output["cancerCellProperties"]["avgPSynthesis"].append(0)
        else:
            hifRates = [a.currentHifRate for a in cancerCells]
            model.output["cancerCellProperties"]["avgHif"].append(
                np.mean(hifRates))

            vegfRates = [a.currentVegfSecretionRate for a in cancerCells]
            model.output["cancerCellProperties"]["avgVegf"].append(
                np.mean(vegfRates))

            metabolicRates = [a.currentMetabolicRate for a in cancerCells]
            model.output["cancerCellProperties"]["avgMetabolicRates"].append(
                np.mean(metabolicRates))

            pSynthesis = [a.currentPSynthesis for a in cancerCells]
            model.output["cancerCellProperties"]["avgPSynthesis"].append(
                np.mean(pSynthesis))

            # Percentage of warburg cells
            warburgCells = [c for c in cancerCells if c.warburgSwitch]
            model.output["cancerCellProperties"]["numWarburgCells"].append(
                float(len(warburgCells)) / float(len(cancerCells)))

            if model.current_epoch % self.distributionInterval == 0 or \
                    model.current_epoch == model.epochs - 1:
                n, bins, patches = plt.hist(hifRates,
                                            bins=list(np.arange(0, 17, 1)))
                model.output["cancerCellProperties"][
                    "HIFExpressionRatesDistributions"].append(
                    {"n": n, "bins": bins, "epoch": model.current_epoch})
