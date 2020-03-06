import numpy as np
import time
from fipy import Grid3D, CellVariable, TransientTerm, DiffusionTerm
from panaxea.core.Steppables import Helper


class GlucoseDiffusionHelper(Helper):
    def __init__(self, model, cancerCellName="CancerCell"):
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.glucoseEnvName = model.properties["envNames"]["glucoseEnvName"]
        self.glucoseDiffusionCoeff = model.properties["diffusion"][
            "glucoseDiffusivity"]
        self.dt = model.properties["diffusion"]["dt"]
        self.diffusionSolveIterations = model.properties["diffusion"][
            "diffusionSolveIterations"]
        self.cancerCellName = cancerCellName
        self.baseGlucoseSecretionRate = \
        model.properties["agents"]["endothelialCells"]["glucoseSecretionRate"]
        self.minGlucoseUptakeRate = model.properties["agents"]["cancerCells"][
            "minGlucoseUptakeRate"]
        self.maxGlucoseUptakeRate = model.properties["agents"]["cancerCells"][
            "maxGlucoseUptakeRate"]

    def getSourceSinkGrids_(self, phi, model):
        start = time.time()
        xsize = model.environments[self.glucoseEnvName].xsize
        ysize = model.environments[self.glucoseEnvName].ysize
        zsize = model.environments[self.glucoseEnvName].zsize

        agentGrid = model.environments["agentEnv"].grid
        phiTmp = np.reshape(phi._array, (xsize, ysize, zsize))

        sourceGrid = CellVariable(name="source",
                                  mesh=Grid3D(dx=1, dy=1, dz=1, nx=xsize,
                                              ny=ysize, nz=zsize))
        sinkGrid = CellVariable(name="sink",
                                mesh=Grid3D(dx=1, dy=1, dz=1, nx=xsize,
                                            ny=ysize, nz=zsize))

        for coordinate, agents in agentGrid.items():
            if len(agents) == 0:
                continue

            concentrationAtPos = phiTmp[coordinate[0]][coordinate[1]][
                coordinate[2]]

            # Getting the current sink rate, defined as the sum of the sink
            # rates of all non-dead and non-quiescent
            # cancer cells and healthy cells at this position
            sinkWarburg = [a.glucoseUptakeRate for a in agents if
                           a.__class__.__name__ == self.cancerCellName
                           and not (a.dead or a.quiescent) and a.warburgSwitch]
            sinkRateWarburg = sum(sinkWarburg)

            sinkNonWarburg = [a.glucoseUptakeRate for a in agents if
                              (a.__class__.__name__ == self.cancerCellName
                               and not (
                                              a.dead or a.quiescent) and not
                               a.warburgSwitch) or (
                                          a.__class__.__name__ ==
                                          "HealthyCell" and not a.dead)]

            sinkRateNonWarburg = sum(sinkNonWarburg)

            # Getting the current source rate, defined as the sum of source
            # rates of all Tip and Trunk cells at
            # this position
            sourceRate = sum(
                [a.glucoseSecretionRate for a in agents if
                 a.__class__.__name__ in ("TipCell", "TrunkCell")])

            # A pre-estimate of what the concentration at this position will
            # be. This of course neglects diffusion,
            # but can give an estimate of how we should regulate our sources
            # and sinks
            estimatedConcentration = concentrationAtPos + sourceRate - \
                                     sinkRateWarburg - sinkRateNonWarburg

            # If our estimated concentration is greater than our source
            # rate, this means we really are outputting
            # too much. At most, we want to achieve equilibrium between
            # sources and environment, so we reduce our
            # output rate. Of course, we can't reduce our output rate by
            # more than the output rate itself
            if estimatedConcentration >= self.baseGlucoseSecretionRate:
                sourceRate -= min(sourceRate,
                                  estimatedConcentration -
                                  self.baseGlucoseSecretionRate)

            # If our estimate concentration is below zero, then our sinks
            # should be reduced. We reduce them by the
            # magnitude of the negative value, but of course we can't reduce
            # them beyond the original value.

            if estimatedConcentration < 0:
                numWarburg = len(sinkWarburg) * 1.
                numNonWarburg = len(sinkNonWarburg) * 1.
                totSink = numWarburg + numNonWarburg * 1.

                ratioWarburg = numWarburg / totSink if totSink > 0 else 0.
                ratioNonWarburg = numNonWarburg / totSink if totSink > 0 \
                    else 0.

                # Sink rates cannot be lower than respective minimum glucose
                # uptake rates, otherwise this means
                # we don't have enough glucose and the cell should die. Each
                # sink will be decreased by the proportion
                # of negative estimate concentration for which they are
                # responsible
                sinkRateWarburg -= min(
                    sinkRateWarburg - self.maxGlucoseUptakeRate * ratioWarburg,
                    abs(estimatedConcentration) * ratioWarburg)
                sinkRateNonWarburg -= min(
                    sinkRateNonWarburg - self.minGlucoseUptakeRate *
                    ratioNonWarburg,
                    abs(estimatedConcentration) * ratioNonWarburg)

            i = np.ravel_multi_index(
                [coordinate[0], coordinate[1], coordinate[2]],
                (xsize, ysize, zsize))

            sourceGrid.value[i] = sourceRate
            sinkGrid.value[i] = sinkRateNonWarburg + sinkRateWarburg

        return sourceGrid, sinkGrid

    def solveDiffusion_(self, model):
        glucoseGrid = model.environments[self.glucoseEnvName]
        nx = glucoseGrid.xsize
        ny = glucoseGrid.ysize
        nz = glucoseGrid.zsize

        dx = dy = dz = 1.0

        D = self.glucoseDiffusionCoeff

        mesh = Grid3D(dx=dx, dy=dy, nx=nx, ny=ny, dz=dz, nz=nz)

        phi = CellVariable(name="solutionvariable", mesh=mesh)
        phi.setValue(0.)

        start = time.time()
        for i in range(self.diffusionSolveIterations):
            sourceGrid, sinkGrid = self.getSourceSinkGrids_(phi, model)
            eq = TransientTerm() == DiffusionTerm(
                coeff=D) + sourceGrid - sinkGrid

            eq.solve(var=phi, dt=1)
            eq = TransientTerm() == DiffusionTerm(coeff=D)

            eq.solve(var=phi, dt=self.dt)
        end = time.time()
        print("Solving glucose diffusion took %s seconds" % str(end - start))

        return phi, nx, ny, nz

    def step_prologue(self, model):

        suitableSolution = False
        iteration = 1
        while not suitableSolution:

            if iteration > 2:
                print(
                            "Glucose diffusion still has negative positions "
                            "at epochs %s despite killing all agents at such "
                            "coordinates..." % str(
                        model.current_epoch))
                print(negativePositions)
                model.exit = True
                break

            print("Solving glucose diffusion")
            print("Solving for iteration %s" % str(iteration))
            phi, nx, ny, nz = self.solveDiffusion_(model)
            cs = phi._array

            cs = np.reshape(cs, (nx, ny, nz))

            # Handling negative concentrations by
            negativePositions = []

            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        if cs[x][y][z] < 0:
                            negativePositions.append(((x, y, z), cs[x][y][z]))

            if len(negativePositions) == 0:
                suitableSolution = True
            else:
                # print("Negative positions (Glucose)")
                # print(negativePositions)

                for p in negativePositions:
                    p = p[0]
                    for a in [a for a in model.environments["agentEnv"].grid[
                        (p[0], p[1], p[2])] if
                              a.__class__.__name__ in ["HealthyCell",
                                                       "CancerCell"]]:
                        a.dead = True

                        if a.__class__.__name__ == "CancerCell":
                            a.causeOfDeath = {
                                "cause": "Lack of glucose",
                                "oxygenAtPos": 0,
                                "warburg": a.warburgSwitch
                            }

                iteration = iteration + 1

        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    model.environments[self.glucoseEnvName].grid[(x, y, z)] = \
                    cs[x][y][z]
