import numpy as np
import time
from fipy import Grid3D, CellVariable, TransientTerm, DiffusionTerm
from panaxea.core.Steppables import Helper


class VegfDiffusionHelper(Helper):
    def __init__(self, model, cancerCellName="CancerCell"):
        self.agent_env_name = model.properties["envNames"]["agentEnvName"]
        self.vegf_env_name = model.properties["envNames"]["vegfEnvName"]
        self.vegf_diffusion_coeff = model.properties["diffusion"][
            "vegfDiffusivity"]
        self.dt = model.properties["diffusion"]["dt"]
        self.diffusion_solve_iterations = model.properties["diffusion"][
            "diffusionSolveIterations"]
        self.cancer_cell_name = cancerCellName
        self.max_vegf = model.properties["agents"]["cancerCells"][
            "maxVegfSecretionRate"]

    def __get_source_sink_grids(self, phi, model):
        xsize = model.environments[self.vegf_env_name].xsize
        ysize = model.environments[self.vegf_env_name].ysize
        zsize = model.environments[self.vegf_env_name].zsize

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

        for coordinate, agents in agent_grid.items():
            if len(agents) == 0:
                continue

            concentration_at_pos = phi_tmp[coordinate[0]][coordinate[1]][
                coordinate[2]]

            source_rate = sum([a.current_vegf_secretion_rate for a in agents if
                               (a.__class__.__name__ == "CancerCell"
                                and not (a.quiescent or a.dead))])
            # A pre-estimate of what the concentration at this position will
            # be. This of course neglects diffusion,
            # but can give an estimate of how we should regulate our sources
            # and sinks
            estimated_concentration = concentration_at_pos + source_rate

            # If our estimated concentration is greater than our maximum
            # source rate, this means we really are outputting
            # too much. At most, we want to achieve equilibrium between
            # sources and environment, so we reduce our
            # output rate. Of course, we can't reduce our output rate by
            # more than the output rate itself
            if estimated_concentration >= self.max_vegf:
                source_rate -= min(source_rate,
                                   estimated_concentration - self.max_vegf)

            i = np.ravel_multi_index(
                [coordinate[0], coordinate[1], coordinate[2]],
                (xsize, ysize, zsize))

            source_grid.value[i] = source_rate

        return source_grid

    def __solve_diffusion(self, model):
        vegf_grid = model.environments[self.vegf_env_name]
        nx = vegf_grid.xsize
        ny = vegf_grid.ysize
        nz = vegf_grid.zsize

        dx = dy = dz = 1.0

        D = self.vegf_diffusion_coeff

        mesh = Grid3D(dx=dx, dy=dy, nx=nx, ny=ny, dz=dz, nz=nz)

        phi = CellVariable(name="solutionvariable", mesh=mesh)
        phi.setValue(0.)

        start = time.time()
        for i in range(self.diffusion_solve_iterations):
            source_grid = self.__get_source_sink_grids(phi, model)
            eq = TransientTerm() == DiffusionTerm(coeff=D) + source_grid

            eq.solve(var=phi, dt=1)
            eq = TransientTerm() == DiffusionTerm(coeff=D)

            eq.solve(var=phi, dt=self.dt)
        end = time.time()
        print("Solving VEGF diffusion took %s seconds" % str(end - start))

        return phi, nx, ny, nz

    def step_prologue(self, model):

        suitable_solution = False
        iteration = 1
        negative_positions = []
        while not suitable_solution:

            if iteration > 2:
                print(
                        "Vegf diffusion still has negative positions at "
                        "epochs %s despite killing all agents at such "
                        "coordinates..." % str(model.current_epoch))
                print(negative_positions)
                model.exit = True
                break

            print("Solving vegf diffusion")
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
                            a.cause_of_death = {
                                "cause": "Lack of oxygen",
                                "oxygenAtPos": 0,
                                "warburg": a.warburg_switch
                            }

                iteration = iteration + 1

        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    model.environments[self.vegf_env_name].grid[(x, y, z)] = \
                        cs[x][y][z]
