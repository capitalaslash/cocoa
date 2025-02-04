#! /usr/bin/env python3

import numpy as np

import cocoa  # type: ignore


def setup(p: cocoa.ProblemFD1D):
    p.name = "fd1d_custom"

    # mesh
    p.start = 0.0  # starting point
    p.n = 101  # nb. of points
    p.h = 1.0 / (p.n - 1)

    # coupling
    p.coupling_type = cocoa.COUPLING_TYPE.medcoupling

    # fields
    p.var_name = "T"
    p.u = np.ones(shape=(p.n,)) * 1.0
    # p.u = pc.VectorFD(p.n)
    p.u_old = p.u.copy()

    # eqn params
    alpha = 1.0
    q_init = 1.0
    p.q = np.ones(shape=(p.n,)) * q_init

    # bcs
    p.bcs = cocoa.FDBCList1D(
        cocoa.FDBC(
            side=cocoa.FD_BC_SIDE.left,
            type=cocoa.FD_BC_TYPE.dirichlet,
            value=0.0,
        ),
        cocoa.FDBC(
            side=cocoa.FD_BC_SIDE.right,
            type=cocoa.FD_BC_TYPE.neumann,
            value=0.5 * q_init / alpha,
        ),
    )

    # time
    p.time = 0.0
    p.final_time = 10.0
    p.dt = 1.0

    # la
    p.m.init(p.n)
    p.rhs = np.zeros(shape=(p.n,))

    # assembly
    def assembly(p: cocoa.ProblemFD1D):
        print("custom assembly")

        rhs = np.zeros(shape=(p.n,))

        for k in range(p.n):
            # diagonal
            p.m.add(k, k, 1.0 + 2.0 * alpha * p.dt / (p.h * p.h))

            # left
            v_left = -alpha * p.dt / (p.h * p.h)
            if k > 0:
                p.m.add(k, k - 1, v_left)
            else:
                # bc: u_-1 = u_1 + 2 * h * bc_value
                p.m.add(k, k + 1, v_left)
                p.bc_start.ghost_values = [v_left]

            # right
            v_right = -alpha * p.dt / (p.h * p.h)
            if k < p.n - 1:
                p.m.add(k, k + 1, v_right)
            else:
                # bc: u_n = u_n-2 + 2 * h * bc_value
                p.m.add(k, k - 1, v_right)
                p.bc_end.ghost_values = [v_right]

            # rhs
            rhs[k] = p.u_old[k] + p.q[k] * p.dt

        p.m.close()
        p.rhs = rhs

    # p.eqn_type = cocoa.EQN_TYPE.heat
    p.set_custom_assembly(assembly)

    print(f"assemblies: {p.assemblies}")

    # io
    p.setup_io("output_custom")

if __name__ == "__main__":
    p = cocoa.ProblemFD1D()
    # fd1d_heat implementation without config file
    # setup(p)
    p.setup(config_file="fd1d_heat.dat")

    p.print()

    while p.run():
        p.advance()
        num_iters = p.solve()
        p.print()

        if num_iters == 0:
            break

    print(f"solution at final point: {p.u[-1]:.16e}")
    assert np.fabs(p.u[-1] - 4.6885564243136019e-06) < 1e-12
