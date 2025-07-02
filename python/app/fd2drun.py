import numpy as np

import cocoa  # type: ignore


def setup(p: cocoa.ProblemFD2D):
    p.name = "fd2d_custom"
    p.debug = False

    # mesh
    p.mesh = cocoa.MeshFD2D((0.0, 0.0), (1.0, 1.0), (10, 10))

    # fields
    p.n_vars = 1
    size = p.mesh.n_pts * p.n_vars
    p.var_names = ["u"]
    u_init = 0.0
    p.u.data = np.ones(shape=(size,)) * u_init
    p.u_old.data = p.u.data.copy()
    q_init = 1.0
    p.q.data = np.ones(shape=(size,)) * q_init
    # p.c[0].data = np.ones((size,)) * 0.0
    # p.c[1].data = np.ones((size,)) * 0.0

    # bcs
    p.bcs = [
        cocoa.FDBCList2D(
            bottom=cocoa.FDBC(
                cocoa.FD_BC_SIDE.bottom,
                cocoa.FD_BC_TYPE.dirichlet,
                cocoa.VectorFD(np.ones((p.mesh.n[0],)) * 0.0),
            ),
            right=cocoa.FDBC(
                cocoa.FD_BC_SIDE.right,
                cocoa.FD_BC_TYPE.neumann,
                cocoa.VectorFD(np.ones((p.mesh.n[1],)) * 0.0),
            ),
            top=cocoa.FDBC(
                cocoa.FD_BC_SIDE.top,
                cocoa.FD_BC_TYPE.neumann,
                cocoa.VectorFD(np.ones((p.mesh.n[0],)) * 5.0),
            ),
            left=cocoa.FDBC(
                cocoa.FD_BC_SIDE.left,
                cocoa.FD_BC_TYPE.neumann,
                cocoa.VectorFD(np.ones((p.mesh.n[1],)) * 0.0),
            ),
        ),
    ]

    # time
    p.time = 0.0
    p.final_time = 40.0
    p.dt = 1.0

    # linear algebra
    p.solver_type = cocoa.FD_SOLVER_TYPE.gauss_seidel
    p.m.init(size, 5)
    p.rhs.resize(size)

    # assembly
    def assembly(p: cocoa.ProblemFD2D):
        nx, ny = p.mesh.n
        alpha = 0.1

        for j in range(ny):
            for i in range(nx):
                id_m = i + j * p.mesh.n[0]
                # x, y = p.mesh.pt([i, j])
                hx, hy = p.mesh.h

                # diagonal
                p.m.add(id_m, id_m, 1.0 + 2.0 * alpha * p.dt * (1 / hx**2 + 1 / hy**2))

                # bottom
                value_bottom = -alpha * p.dt / hy**2
                if j > 0:
                    p.m.add(id_m, id_m - nx, value_bottom)
                else:
                    p.m.add(id_m, id_m + nx, value_bottom)
                    p.bcs[0].bottom.ghost_values.set(i, value_bottom)

                # right
                value_right = -alpha * p.dt / hx**2
                if i < nx - 1:
                    p.m.add(id_m, id_m + 1, value_right)
                else:
                    p.m.add(id_m, id_m - 1, value_right)
                    p.bcs[0].right.ghost_values.set(j, value_right)

                # top
                value_top = -alpha * p.dt / hy**2
                if j < ny - 1:
                    p.m.add(id_m, id_m + nx, value_top)
                else:
                    p.m.add(id_m, id_m - nx, value_top)
                    p.bcs[0].top.ghost_values.set(i, value_top)

                # left
                value_left = -alpha * p.dt / hx**2
                if i > 0:
                    p.m.add(id_m, id_m - 1, value_left)
                else:
                    p.m.add(id_m, id_m + 1, value_left)
                    p.bcs[0].left.ghost_values.set(j, value_left)

                # rhs
                p.rhs.add(id_m, p.u_old[id_m] + p.q[id_m])

        p.m.close()

    p.set_custom_assembly(assembly)

    # p.eqn_type = cocoa.EQN_TYPE.heat
    # p.params.set_scalar("alpha", 0.1)

    p.setup_io("output_fd2d_custom")


def main():
    p = cocoa.ProblemFD2D()

    # p.setup(config_file="fd2d_heat.dat")
    setup(p)

    p.print()

    while p.run():
        p.advance()
        num_iters = p.solve()
        p.print()

        if num_iters == 0:
            break


if __name__ == "__main__":
    main()
