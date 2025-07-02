#! /usr/bin/env python3

import numpy as np

import cocoa  # type: ignore


def setup(p: cocoa.ProblemFD1D):
    p.name = "fd1d2v"
    p.debug = False

    # mesh
    p.mesh = cocoa.MeshFD1D(start=[0.0], end=[1.0], n=[10])

    # fields
    p.n_vars = 2
    p.var_names = ["u", "v"]
    size = p.n_vars * p.mesh.n_pts
    p.u.data = np.ones(shape=(size,)) * 0.0
    p.u_old.data = p.u.data.copy()

    # eqn params
    alpha = [0.1, 0.01]
    q_init = [0.0, 0.0]
    q = np.ones(shape=(size,))
    for v in range(p.n_vars):
        start = 0 + v * p.mesh.n_pts
        end = (v + 1) * p.mesh.n_pts
        q[start:end] *= q_init[v]

    # bcs
    p.bcs = [
        # u
        cocoa.FDBCList1D(
            left=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.left,
                type=cocoa.FD_BC_TYPE.dirichlet,
                value=1.0,
            ),
            right=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.right,
                type=cocoa.FD_BC_TYPE.neumann,
                value=0.0,
            ),
        ),
        # v
        cocoa.FDBCList1D(
            left=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.left,
                type=cocoa.FD_BC_TYPE.dirichlet,
                value=1.0,
            ),
            right=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.right,
                type=cocoa.FD_BC_TYPE.neumann,
                value=0.0,
            ),
        ),
    ]

    # time
    p.time = 0.0
    p.final_time = 20.0
    p.dt = 1.0

    # la
    p.max_iters = 10_000
    p.m.init(size, 4)
    p.rhs.resize(size)

    # assembly
    def assembly(p: cocoa.ProblemFD1D):
        print("custom assembly")

        h = p.mesh.h[0]

        for v in range(p.n_vars):
            for k in range(p.mesh.n_pts):
                id_m = k + v * p.mesh.n_pts

                # diagonal
                p.m.add(id_m, id_m, 1.0 + 2.0 * alpha[v] * p.dt / (h * h))

                # left
                value_left = -alpha[v] * p.dt / (h * h)
                if k > 0:
                    p.m.add(id_m, id_m - 1, value_left)
                else:
                    p.m.add(id_m, id_m + 1, value_left)
                    p.bcs[v].left.ghost_values.set(0, value_left)

                # right
                value_right = -alpha[v] * p.dt / (h * h)
                if k < p.mesh.n_pts - 1:
                    p.m.add(id_m, id_m + 1, value_right)
                else:
                    p.m.add(id_m, id_m - 1, value_right)
                    p.bcs[v].right.ghost_values.set(0, value_right)

                # coupling
                id_o = k + ((v + 1) % p.n_vars) * p.mesh.n_pts
                p.m.add(id_m, id_o, 1.0e-1 * p.dt)

                # rhs
                p.rhs.set(id_m, p.u_old[id_m] + q[id_m] * p.dt)

        p.m.close()

    # p.eqn_type = cocoa.EQN_TYPE.heat
    p.set_custom_assembly(assembly)

    # io
    p.setup_io("output_fd1d2v")

    # p.debug = True


if __name__ == "__main__":
    p = cocoa.ProblemFD1D()

    setup(p)

    p.print()

    while p.run():
        p.advance()
        num_iters = p.solve()
        p.print()

        if num_iters == 0:
            break

    # print(f"solution at final point: {p.u.data[-1]:.16e}")
    # assert np.fabs(p.u.data[-1] - TEMP_END) < 10 * SOLVER_TOL
