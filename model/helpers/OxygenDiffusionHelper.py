import numpy as np
import time
from fipy import Grid3D, CellVariable, TransientTerm, DiffusionTerm
from panaxea.core.Steppables import Helper


class OxygenDiffusionHelper(Helper):
    def __init__(self, model, cancerCellName="CancerCell"):
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.oxygen_env_name = model.properties["envNames"]["oxygenEnvName"]
        self.oxygen_diffusion_coeff = model.properties["diffusion"][
            "oxygenDiffusivity"]
        self.dt = model.properties["diffusion"]["dt"]
        self.diffusion_solve_iterations = model.properties["diffusion"][
            "diffusionSolveIterations"]
        self.cancer_cell_name = cancerCellName
        self.base_oxygen_emission_rate = \
            model.properties["agents"]["endothelialCells"][
                "baseOxygenEmissionRate"]
        self.minimun_oxygen_cancer_cells = \
            model.properties["agents"]["cancerCells"][
                "minimumOxygenConcentration"]

    def __get_source_sink_grids(self, phi, model):
        xsize = model.environments[self.oxygen_env_name].xsize
        ysize = model.environments[self.oxygen_env_name].ysize
        zsize = model.environments[self.oxygen_env_name].zsize

        agent_grid = model.environments["agentEnv"].grid
        phi_tmp = np.reshape(phi._array, (xsize, ysize, zsize))

        source_grid = CellVariable(
            name="source",
            mesh=Grid3D(
                dx=1,
                dy=1,
                dz=1,
                nx=xsize,
                ny=ysize,
                nz=zsize)
        )

        sink_grid = CellVariable(
            name="sink",
            mesh=Grid3D(
                dx=1,
                dy=1,
                dz=1,
                nx=xsize,
                ny=ysize,
                nz=zsize)
        )

        for coordinate, agents in agent_grid.items():
            if len(agents) == 0:
                continue

            concentration_at_pos = phi_tmp[coordinate[0]][coordinate[1]][
                coordinate[2]]

            # Getting the current sink rate, defined as the sum of the sink
            # rates of all non-dead and non-quiescent
            # cancer cells and healthy cells at this position
            sink_rate = sum([a.current_metabolic_rate for a in agents if
                             (a.__class__.__name__ == "CancerCell"
                              and not (a.quiescent or a.dead)) or (
                                     a.__class__.__name__ ==
                                     "HealthyCell" and not a.dead)])

            # Getting the current source rate, defined as the sum of source
            # rates of all Tip and Trunk cells at
            # this position
            source_rate = sum(
                [a.oxygen_emission_rate for a in agents if
                 a.__class__.__name__ in ("TipCell", "TrunkCell")])

            # A pre-estimate of what the concentration at this position will
            # be. This of course neglects diffusion,
            # but can give an estimate of how we should regulate our sources
            # and sinks

            rate_diff = source_rate - sink_rate

            estimated_concentration = concentration_at_pos + rate_diff

            # If our estimated concentration is greater than our source
            # rate, this means we really are outputting
            # too much. At most, we want to achieve equilibrium between
            # sources and environment, so we reduce our
            # output rate. Of course, we can't reduce our output rate by
            # more than the output rate itself
            if estimated_concentration >= self.base_oxygen_emission_rate:
                source_rate -= min(source_rate,
                                   estimated_concentration -
                                   self.base_oxygen_emission_rate)

            # If our estimate concentration is below zero, then our sinks
            # should be reduced. We reduce them by the
            # magnitude of the negative value, but of course we can't reduce
            # them beyond the original value.

            if estimated_concentration < 0:
                sink_rate -= min(sink_rate - self.minimun_oxygen_cancer_cells,
                                 abs(estimated_concentration))

            i = np.ravel_multi_index(
                [coordinate[0], coordinate[1], coordinate[2]],
                (xsize, ysize, zsize))

            source_grid.value[i] = source_rate
            sink_grid.value[i] = sink_rate

        return source_grid, sink_grid

    def __solve_diffusion(self, model):
        oxygen_grid = model.environments[self.oxygen_env_name]
        nx = oxygen_grid.xsize
        ny = oxygen_grid.ysize
        nz = oxygen_grid.zsize

        dx = dy = dz = 1.0

        D = self.oxygen_diffusion_coeff

        mesh = Grid3D(dx=dx, dy=dy, nx=nx, ny=ny, dz=dz, nz=nz)

        phi = CellVariable(name="solutionvariable", mesh=mesh)
        phi.setValue(0.)

        start = time.time()
        for _ in range(self.diffusion_solve_iterations):
            source_grid, sink_grid = self.__get_source_sink_grids(phi, model)
            eq = TransientTerm() == DiffusionTerm(
                coeff=D) + source_grid - sink_grid

            eq.solve(var=phi, dt=1)
            eq = TransientTerm() == DiffusionTerm(coeff=D)

            eq.solve(var=phi, dt=self.dt)
        end = time.time()
        print("Solving oxygen diffusion took %s seconds" % str(end - start))

        return phi, nx, ny, nz

    def step_prologue(self, model):
        suitable_solution = False
        iteration = 1
        negative_positions = []

        while not suitable_solution:

            if iteration > 2:
                print(
                        "Oxygen diffusion still has negative positions "
                        "at epochs %s despite killing all agents at such "
                        "coordinates..." % str(model.current_epoch))
                print(negative_positions)
                model.exit = True
                break

            print("Solving oxygen diffusion")
            print("Solving for iteration %s" % str(iteration))
            phi, nx, ny, nz = self.__solve_diffusion(model)
            cs = phi._array

            cs = np.reshape(cs, (nx, ny, nz))

            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        if cs[x][y][z] < 0:
                            negative_positions.append(((x, y, z), cs[x][y][z]))

            if len(negative_positions) == 0:
                suitable_solution = True
            else:
                # print("Negative positions (Oxygen)")
                # print(negativePositions)

                for p in negative_positions:
                    p = p[0]
                    for a in [a for a in model.environments["agentEnv"].grid[
                        (p[0], p[1], p[2])] if
                              a.__class__.__name__ in ["HealthyCell",
                                                       "CancerCell"]]:
                        a.dead = True

                        if a.__class__.__name__ == "CancerCell":
                            a.causeOfDeath = {
                                "cause": "Lack of oxygen",
                                "oxygenAtPos": 0,
                                "warburg": a.warburg_switch
                            }

                iteration = iteration + 1

        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    model.environments[self.oxygen_env_name].grid[(x, y,
                                                                   z)] = \
                        cs[x][y][z]
