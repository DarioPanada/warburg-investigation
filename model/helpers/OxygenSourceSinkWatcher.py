from panaxea.core.Steppables import Helper


class OxygenSourceSinkWatcher(object, Helper):
    '''
    Will store the source/sink coordinates for oxygen as part of the model
    properties

    :param model - The main model object
    :param oxygenDiffusionHelper - A reference to the oxygen diffusion helper
    :param interval - (Optional - defaults to 10), will store every
    'interval' epochs
    '''

    def __init__(self, model, oxygenDiffusionHelper, interval=10):
        model.output["oxygenGrids"] = []
        self.oxygen_diffusion_helper = oxygenDiffusionHelper
        self.interval = interval

    def step_epilogue(self, model):
        # Diffusion runs in the prologue, so we capture in the epilogue.
        if model.current_epoch % self.interval == 0:
            source_coords = self.oxygen_diffusion_helper.sourceCoords
            sink_coords = self.oxygen_diffusion_helper.sinkCoords
            model.output["oxygenGrids"].append({
                "epoch": model.current_epoch,
                "sourceCoords": source_coords,
                "sinkCoords": sink_coords
            })
