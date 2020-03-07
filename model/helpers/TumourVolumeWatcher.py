from cmath import sqrt

from panaxea.core.Steppables import Helper

'''
At each time-step, calculates an estimate of the tumour volume by looking at 
the euclidean distance of the two farthest
cancer cells.
'''


class TumourVolumeWatcher(Helper, object):
    def __init__(self, model, cancerCellClassName="CancerCell"):
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.cancerCellClassName = cancerCellClassName

        model.output["maxDistances"] = []

    def step_epilogue(self, model):
        cancerCellsCoords = [a.environment_positions[self.agentEnvName] for a in
                             model.schedule.agents if
                             a.__class__.__name__ == self.cancerCellClassName]
        scoredCoords = [(sum(c), c) for c in cancerCellsCoords]

        def f(c):
            return c[0]

        try:
            maxCoord = max(scoredCoords, key=f)[1]
            minCoord = min(scoredCoords, key=f)[1]

            dist = sqrt((maxCoord[0] - minCoord[0]) ** 2 + (
                        maxCoord[1] - minCoord[1]) ** 2 + (
                                    maxCoord[2] - minCoord[2]) ** 2)
        except:
            dist = 0

        model.output["maxDistances"].append(dist)
