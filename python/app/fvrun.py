#! /usr/bin/env python3

import numpy as np

import cocoa  # type: ignore


def setup(p: cocoa.ProblemFD1D):
    p.name = "fv1d"
    p.debug = False

    # mesh
    p.mesh = cocoa.MeshFD1D(start=[0.0], end=[1.0], n=[10])

    # fields
    p.n_vars = 1
    p.var_names = ["u"]
    p.u.data = np.ones(shape=(p.mesh.n_pts,)) * 0.0
    p.u_old.data = p.u.data.copy()

    # io
    p.setup_io("output_fv1d")

    # eqn params
    alpha = np.ones(shape=(p.mesh.n_pts,))
    alpha[:6] *= 1.0
    alpha[6:] *= 0.1
    q = np.ones(shape=(p.mesh.n_pts,)) * 1.0

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
                value=0.0,
            ),
        )
    ]

    # time
    p.time = 0.0
    p.final_time = 20.0
    p.dt = 1.0

    # la
    p.solver_type = cocoa.FD_SOLVER_TYPE.cg
    p.max_iters = 10_000
    p.tol = 1e-6
    p.m.init(p.mesh.n_pts, 3)
    p.rhs.resize(p.mesh.n_pts)

    # assembly
    def assembly(p: cocoa.ProblemFD1D):
        print("custom assembly")

        n = p.mesh.n_pts
        h = p.mesh.h[0]
        ih2 = 1 / (h * h)
        idt = 1 / p.dt

        for k in range(0, n):

            # diagonal
            # p.m.add(k, k, idt + 2.0 * alpha * ih2)
            # p.m.add(k, k, idt + 2 * alpha[k] * ih2)
            if 0 < k < n - 1:
                p.m.add(
                    k, k, idt + 0.5 * (alpha[k - 1] + 2 * alpha[k] + alpha[k + 1]) * ih2
                )
            else:
                p.m.add(k, k, idt + 2 * alpha[k] * ih2)

            # left
            # value_left = -alpha * ih2
            value_left = -alpha[k] * ih2
            if 0 < k < n - 1:
                p.m.add(k, k - 1, -0.5 * (alpha[k - 1] + alpha[k]) * ih2)
            elif k == n - 1:
                p.m.add(k, k - 1, value_left)
            else:  # k == 0
                p.m.add(k, k + 1, value_left)
                p.bcs[0].left.ghost_values.set(0, value_left)

            # right
            # value_right = -alpha * ih2
            value_right = -alpha[k] * ih2
            if 0 < k < n - 1:
                p.m.add(k, k + 1, -0.5 * (alpha[k + 1] + alpha[k]) * ih2)
            elif k == 0:
                p.m.add(k, k + 1, value_right)
            else:  # k == n - 1
                p.m.add(k, k - 1, value_right)
                p.bcs[0].right.ghost_values.set(0, value_right)

            # rhs
            p.rhs.set(k, p.u_old[k] * idt + q[k])

        p.m.close()

    p.set_custom_assembly(assembly)


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
