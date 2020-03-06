from collections import defaultdict

import time
from core.Steppables import Helper
from fipy import Grid3D, CellVariable, TransientTerm, DiffusionTerm
import numpy as np


class DrugDiffusionHelper(object, Helper):
    def __init__(self, model, cancerCellName="CancerCell"):
        self.agentEnvName = model.properties["envNames"]["agentEnvName"]
        self.drugEnvName = model.properties["envNames"]["drugEnvName"]
        self.drugDiffusionCoeff = model.properties["drug"]["diffusivity"]
        self.drugName = model.properties["drug"]["name"]
        self.dt = model.properties["diffusion"]["dt"]
        self.diffusionSolveIterations = model.properties["diffusion"]["diffusionSolveIterations"]
        self.cancerCellName = cancerCellName

    def stepPrologue(self, model):
        # Coordinates of sinks (cancer cells)
        sinkAgents = [a for a in model.schedule.agents if a.__class__.__name__ == self.cancerCellName and not a.dead and not a.quiescent]
        sinkCoords = defaultdict(int)

        for a in sinkAgents:
            sinkCoords[a.environmentPositions[self.agentEnvName]] += a.drugUptakeRate

        # Coordinates of sources (endothelial cells)
        sourceAgents = [a for a in model.schedule.agents if
                        a.__class__.__name__ == "TipCell" or a.__class__.__name__ == "TrunkCell"]
        sourceCoords = defaultdict(int)

        for a in sourceAgents:
            sourceCoords[a.environmentPositions[self.agentEnvName]] += a.drugSecretionRate

        drugGrid = model.environments[self.drugEnvName]
        nx = drugGrid.xsize
        ny = drugGrid.ysize
        nz = drugGrid.zsize

        dx = dy = dz = 1.0

        D = self.drugDiffusionCoeff

        mesh = Grid3D(dx=dx, dy=dy, nx=nx, ny=ny, dz=dz, nz=nz)

        phi = CellVariable(name="solutionvariable", mesh=mesh, value=0.)
        phi.setValue(0.)

        x, y, z = mesh.cellCenters

        sourceGrid = self.setupSourceGrid_(sourceCoords, mesh, x, y, z)
        sinkGrid = self.setupSinkGrid_(sinkCoords, mesh, x, y, z)

        eq = TransientTerm() == DiffusionTerm(coeff=D) + sourceGrid - sinkGrid

        # solver = DefaultSolver(iterations=200)

        start = time.time()
        for _ in self.diffusionSolveIterations:
            eq.solve(var=phi, dt=self.dt)
        end = time.time()
        print("Solving drug diffusion took %s seconds" % str(end - start))

        cs = phi._array
        cs = np.reshape(cs, (nx, ny, nz))

        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    model.environments[self.drugEnvName].grid[(x, y, z)] = cs[x][y][z]

    def setupSourceGrid_(self, sourceCoords, mesh, x, y, z, patchSize=1):
        sourceGrid = CellVariable(name="source", mesh=mesh, value=0)
        sourceGrid.setValue(0.)
        for pos, v in sourceCoords.iteritems():
            sourceGrid.setValue(v,
                                where=(z > pos[0] - patchSize) & (z < pos[0] + patchSize) & (y > pos[1] - patchSize) & (
                                y < pos[1] + patchSize) & (x > pos[2] - patchSize) & (x < pos[2] + patchSize))

        return sourceGrid

    def setupSinkGrid_(self, sinkCoords, mesh, x, y, z, patchSize=1):
        sinkGrid = CellVariable(name="source", mesh=mesh, value=0)
        sinkGrid.setValue(0.)
        for pos, v in sinkCoords.iteritems():
            sinkGrid.setValue(v,
                              where=(z > pos[0] - patchSize) & (z < pos[0] + patchSize) & (y > pos[1] - patchSize) & (
                              y < pos[1] + patchSize) & (x > pos[2] - patchSize) & (x < pos[2] + patchSize))

        return sinkGrid
