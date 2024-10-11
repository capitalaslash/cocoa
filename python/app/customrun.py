#! /usr/bin/env python3

import numpy as np

import pycocoa as pc

if __name__ == "__main__":
    p = pc.ProblemFD1D()

    # replace setup from file with local setup
    # p.setup(config_file="fd1d_heat.dat")

    p.name = "fd1d_heat"

    # mesh
    p.start = 0.0  # starting point
    p.n = 6  # nb. of points
    p.h = 1.0 / (p.n - 1)

    # coupling
    p.couplingType = pc.COUPLING_TYPE.medcoupling  # pc.COUPLING_TYPE.simple

    # fields
    p.u = [1.0] * p.n  # pc.VectorFD(p.n)  # np.ndarray(shape=[1, p.n])
    # tmp = pc.VectorFD(p.n)
    # print(type(tmp))
    # print(tmp)
    # print(type(p.u))
    # p.u = pc.VectorFD(p.n)
    p.uOld = [0.0] * p.n

    # eqn params
    p.q = [20.0] * p.n
    p.alpha = 2.0

    # time
    p.time = 0.0
    p.finalTime = 4.0
    p.dt = 0.4

    # la
    p.m.init(p.n)
    p.rhs = [0.0] * p.n

    # assembly
    p.eqnType = pc.EQN_TYPE.heat

    # bcs
    p.bcStart = pc.FDBC(pc.FD_BC_TYPE.dirichlet, 1.0)
    p.bcEnd = pc.FDBC(pc.FD_BC_TYPE.neumann, 0.0)

    # io
    p.setupIO("output_custom")

    p.print()

    while p.run():
        p.advance()
        p.solve()
        p.print()
