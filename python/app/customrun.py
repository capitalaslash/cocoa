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
    p.q = np.ones(shape=(p.n,)) * 1.0
    p.alpha = 1.0

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

        # eqn parameters
        alpha = 0.1

        rhs = np.zeros(shape=(p.n,))

        for k in range(p.n):
            kLeft = k - 1 if k != 0 else k + 1
            kRight = k + 1 if k != p.n - 1 else k - 1

            p.m.add(k, k, 1.0 + 2.0 * alpha * p.dt / (p.h * p.h))
            p.m.add(k, kLeft, -alpha * p.dt / (p.h * p.h))
            p.m.add(k, kRight, -alpha * p.dt / (p.h * p.h))
            rhs[k] = p.uOld[k] + p.q[k]

        p.m.close()
        p.rhs = rhs

    p.setAssemblyCustom(assembly)

    # bcs
    p.bcStart = pc.FDBC(pc.FD_BC_TYPE.dirichlet, [1.0])
    p.bcEnd = pc.FDBC(pc.FD_BC_TYPE.neumann, [0.0])

    # io
    p.setupIO("output_custom")

if __name__ == "__main__":
    # fd1d_heat implementation without config file
    p = pc.ProblemFD1D()
    setup(p)

    p.print()

    while p.run():
        p.advance()
        p.solve()
        p.print()
