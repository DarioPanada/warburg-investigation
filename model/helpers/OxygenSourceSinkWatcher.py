from core.Steppables import Helper

'''
Will store the source/sink coordinates for oxygen as part of the model properties
'''


class OxygenSourceSinkWatcher(object, Helper):
    '''
    :param model - The main model object
    :param oxygenDiffusionHelper - A reference to the oxygen diffusion helper
    :param interval - (Optional - defaults to 10), will store every 'interval' epochs
    '''

    def __init__(self, model, oxygenDiffusionHelper, interval=10):
        model.output["oxygenGrids"] = []
        self.oxygenDiffusionHelper = oxygenDiffusionHelper
        self.interval = interval

    '''
    Diffusion runs in the prologue, so we capture in the epilogue.
    '''

    def stepEpilogue(self, model):
        if model.currentEpoch % self.interval == 0:
            sourceCoords = self.oxygenDiffusionHelper.sourceCoords
            sinkCoords = self.oxygenDiffusionHelper.sinkCoords
            model.output["oxygenGrids"].append({
                "epoch": model.currentEpoch,
                "sourceCoords": sourceCoords,
                "sinkCoords": sinkCoords
            })
