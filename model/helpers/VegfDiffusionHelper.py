from collections import defaultdict

import time
from core.Steppables import Helper
from fipy import Grid3D, CellVariable, TransientTerm, DiffusionTerm
import numpy as np


class VegfDiffusionHelper(object, Helper):
    def __init__(self, model, cancerCellName="CancerCell"):
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.vegfEnvName = model.properties["envNames"]["vegfEnvName"]
        self.vegfDiffusionCoeff = model.properties["diffusion"]["vegfDiffusivity"]
        self.dt = model.properties["diffusion"]["dt"]
        self.diffusionSolveIterations = model.properties["diffusion"]["diffusionSolveIterations"]
        self.cancerCellName = cancerCellName
        self.maxVegf = model.properties["agents"]["cancerCells"]["maxVegfSecretionRate"]


    def getSourceSinkGrids_(self, phi, model):
        start = time.time()
        xsize = model.environments[self.vegfEnvName].xsize
        ysize = model.environments[self.vegfEnvName].ysize
        zsize = model.environments[self.vegfEnvName].zsize

        agentGrid = model.environments["agentEnv"].grid
        phiTmp = np.reshape(phi._array, (xsize, ysize, zsize))

        sourceGrid = CellVariable(name="source", mesh=Grid3D(dx=1, dy=1, dz=1, nx=xsize, ny=ysize, nz=zsize))

        for coordinate, agents in agentGrid.items():
            if len(agents) == 0:
                continue

            concentrationAtPos = phiTmp[coordinate[0]][coordinate[1]][coordinate[2]]

            sourceRate = sum([a.currentVegfSecretionRate for a in agents if (a.__class__.__name__ == "CancerCell"
                                                                       and not (a.quiescent or a.dead))])
            # A pre-estimate of what the concentration at this position will be. This of course neglects diffusion,
            # but can give an estimate of how we should regulate our sources and sinks
            estimatedConcentration = concentrationAtPos + sourceRate

            # If our estimated concentration is greater than our maximum source rate, this means we really are outputting
            # too much. At most, we want to achieve equilibrium between sources and environment, so we reduce our
            # output rate. Of course, we can't reduce our output rate by more than the output rate itself
            if estimatedConcentration >= self.maxVegf:
                sourceRate -= min(sourceRate, estimatedConcentration - self.maxVegf)

            i = np.ravel_multi_index([coordinate[0], coordinate[1], coordinate[2]], (xsize, ysize, zsize))

            sourceGrid.value[i] = sourceRate

        return sourceGrid

    def solveDiffusion_(self, model):
        vegfGrid = model.environments[self.vegfEnvName]
        nx = vegfGrid.xsize
        ny = vegfGrid.ysize
        nz = vegfGrid.zsize

        dx = dy = dz = 1.0

        D = self.vegfDiffusionCoeff

        mesh = Grid3D(dx=dx, dy=dy, nx=nx, ny=ny, dz=dz, nz=nz)

        phi = CellVariable(name="solutionvariable", mesh=mesh)
        phi.setValue(0.)

        start = time.time()
        for i in range(self.diffusionSolveIterations):
            sourceGrid = self.getSourceSinkGrids_(phi, model)
            eq = TransientTerm() == DiffusionTerm(coeff=D) + sourceGrid

            eq.solve(var=phi, dt=1)
            eq = TransientTerm() == DiffusionTerm(coeff=D)

            eq.solve(var=phi, dt=self.dt)
        end = time.time()
        print("Solving VEGF diffusion took %s seconds" % str(end - start))

        return phi, nx, ny, nz

    def stepPrologue(self, model):

        suitableSolution = False
        iteration = 1
        while not suitableSolution:

            if iteration > 2:
                print("Vegf diffusion still has negative positions at epochs %s despite killing all agents at such coordinates..." % str(model.currentEpoch))
                print(negativePositions)
                model.exit = True
                break

            print("Solving vegf diffusion")
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
                #print("Negative positions (Oxygen)")
                #print(negativePositions)

                for p in negativePositions:
                    p = p[0]
                    for a in [a for a in model.environments["agentEnv"].grid[(p[0], p[1], p[2])] if
                              a.__class__.__name__ in ["HealthyCell", "CancerCell"]]:
                        a.dead = True

                        if a.__class__.__name__ == "CancerCell":
                            a.causeOfDeath = {
                                "cause": "Lack of oxygen",
                                "oxygenAtPos": 0,
                                "warburg": a.warburgSwitch
                            }

                iteration = iteration + 1

        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    model.environments[self.vegfEnvName].grid[(x, y, z)] = cs[x][y][z]

