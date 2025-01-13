#! /usr/bin/env python3

import numpy as np

import pycocoa as pc  # type: ignore


def setup(p: pc.ProblemFD1D):
    p.name = "fd1d_custom"

    # mesh
    p.start = 0.0  # starting point
    p.n = 101  # nb. of points
    p.h = 1.0 / (p.n - 1)

    # coupling
    p.couplingType = pc.COUPLING_TYPE.medcoupling

    # fields
    p.varName = "T"
    p.u = np.ones(shape=(p.n,)) * 1.0
    # p.u = pc.VectorFD(p.n)
    p.uOld = p.u.copy()

    # eqn params
    alpha = 1.0
    q_init = 1.0
    p.q = np.ones(shape=(p.n,)) * q_init

    # bcs
    p.bc_start = pc.FDBC(pc.FD_BC_TYPE.dirichlet, [0.0])
    p.bc_end = pc.FDBC(pc.FD_BC_TYPE.neumann, [0.5 * q_init / alpha])

    # time
    p.time = 0.0
    p.finalTime = 100.0
    p.dt = 0.5

    # la
    p.m.init(p.n)
    p.rhs = np.zeros(shape=(p.n,))

    # assembly
    # p.eqnType = pc.EQN_TYPE.heat

    def assembly(p: pc.ProblemFD1D):
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

    p.setAssemblyCustom(assembly)

    # bcs
    p.bcStart = pc.FDBC(pc.FD_BC_TYPE.dirichlet, [1.0])
    p.bcEnd = pc.FDBC(pc.FD_BC_TYPE.neumann, [0.0])

    # io
    p.setupIO("output_custom")

if __name__ == "__main__":
    p = pc.ProblemFD1D()
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
