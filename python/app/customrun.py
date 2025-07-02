#! /usr/bin/env python3

import numpy as np

import cocoa  # type: ignore


TEMP_END = 2.0
SOLVER_TOL = 1.0e-8


def setup(p: cocoa.ProblemFD1D):
    p.name = "fd1d_custom"
    p.debug = False

    # mesh
    p.mesh = cocoa.MeshFD1D(start=[0.0], end=[1.0], n=[10])

    # fields
    p.n_vars = 1
    p.var_names = ["T"]
    p.u.data = np.ones(shape=(p.mesh.n_pts,)) * 1.0
    p.u_old.data = p.u.data.copy()

    # eqn params
    alpha = 1.0
    q_init = 1.0
    q = np.ones(shape=(p.mesh.n_pts,)) * q_init

    # bcs
    p.bcs = [
        cocoa.FDBCList1D(
            left=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.left,
                type=cocoa.FD_BC_TYPE.dirichlet,
                value=0.0,
            ),
            right=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.right,
                type=cocoa.FD_BC_TYPE.neumann,
                value=TEMP_END - 0.5 * q_init / alpha,
            ),
        )
    ]

    # time
    p.time = 0.0
    p.final_time = 20.0
    p.dt = 1.0

    # la
    p.max_iters = 10_000
    p.tol = SOLVER_TOL
    p.m.init(p.mesh.n_pts, 3)
    p.rhs.resize(p.mesh.n_pts)

    # assembly
    def assembly(p: cocoa.ProblemFD1D):
        print("custom assembly")

        h = p.mesh.h[0]

        for k in range(p.mesh.n_pts):

            # diagonal
            p.m.add(k, k, 1.0 + 2.0 * alpha * p.dt / (h * h))

            # left
            value_left = -alpha * p.dt / (h * h)
            if k > 0:
                p.m.add(k, k - 1, value_left)
            else:
                # bc: u_-1 = u_1 + 2 * h * bc_value
                p.m.add(k, k + 1, value_left)
                p.bcs[0].left.ghost_values.set(0, value_left)

            # right
            value_right = -alpha * p.dt / (h * h)
            if k < p.mesh.n_pts - 1:
                p.m.add(k, k + 1, value_right)
            else:
                # bc: u_n = u_n-2 + 2 * h * bc_value
                p.m.add(k, k - 1, value_right)
                p.bcs[0].right.ghost_values.set(0, value_right)

            # rhs
            p.rhs.set(k, p.u_old[k] + q[k] * p.dt)

        p.m.close()

    # p.eqn_type = cocoa.EQN_TYPE.heat
    p.set_custom_assembly(assembly)

    print(f"assemblies: {p.assemblies}")

    # io
    p.setup_io("output_custom")


if __name__ == "__main__":
    p = cocoa.ProblemFD1D()

    # fd1d_heat implementation without config file
    setup(p)

    p.print()

    while p.run():
        p.advance()
        num_iters = p.solve()
        p.print()

        if num_iters == 0:
            break

    print(f"solution at final point: {p.u.data[-1]:.16e}")
    assert np.fabs(p.u.data[-1] - TEMP_END) < 10 * SOLVER_TOL
