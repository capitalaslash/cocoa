import numpy as np

import cocoa  # type: ignore
from fdutils import Region, Target, TargetFunc

# == USER INPUT ==

N_ELEMS = (20, 20)
ALPHA = 1.0e-1
BETA = 1.0e-3
T_INIT = 0.0
FLUX_EXT = 0.5
DT = 1.0e-1
FINAL_TIME = 100.0
PRINT_STEP = 10
FUNC_TOL = 1e-4

# TEMP_TARGETS = [
#     Target(value=1.0, left=0.4, right=0.8, bottom=0.2, top=0.8),
# ]
TEMP_TARGET = TargetFunc(
    left=0.4,
    right=0.8,
    bottom=0.1,
    top=0.9,
    values=np.array([1.0, 0.7, 0.6, 0.9, 0.9, 0.8, 0.8, 1.1, 1.0]),
)

CONTROL_REGIONS = [
    Region(left=0.1, right=0.3, bottom=0.1, top=0.9),
]

DOMAIN = Region(left=0.0, right=1.0, bottom=0.0, top=1.0)

# == USER INPUT END ==

H = (
    (DOMAIN.right - DOMAIN.left) / N_ELEMS[0],
    (DOMAIN.top - DOMAIN.bottom) / N_ELEMS[1],
)

EPS = 0.5 * min(H[0], H[1])
for c in CONTROL_REGIONS:
    c.eps = EPS
# for t in TEMP_TARGETS:
#     t.eps = EPS
TEMP_TARGET.eps = EPS


def setup_oc(p: cocoa.ProblemFD2D):
    p.name = "fd2d_oc_distributed"
    p.max_iters = 10_000

    # mesh
    p.mesh = cocoa.MeshFD2D(
        start=(DOMAIN.left, DOMAIN.bottom), end=(DOMAIN.right, DOMAIN.top), n=N_ELEMS
    )

    # coupling
    p.coupling_type = cocoa.COUPLING_TYPE.medcoupling

    # fields
    p.n_vars = 2
    size = p.mesh.n_pts * p.n_vars
    p.var_names = ["u", "v"]
    p.u.data = np.ones(shape=(size,)) * T_INIT
    p.u_old.data = p.u.data.copy()

    # parameters
    p.params.set_scalar("alpha", ALPHA)
    p.params.set_scalar("beta", BETA)
    # p.params.set_vector("target", TEMP_TARGETS[0].to_list())
    p.params.set_vector("target", TEMP_TARGET.to_list())
    p.params.set_vector("control", CONTROL_REGIONS[0].to_list())

    # io
    p.setup_io("output_fd2d_custom")
    p.print_step = PRINT_STEP
    # p.clean_output = True

    # discretized target
    target_data = cocoa.VectorFD(p.mesh.n_pts, 0.0)
    for j in range(p.mesh.n[1]):
        for i in range(p.mesh.n[0]):
            pt_id = i + p.mesh.n[0] * j
            pt = p.mesh.pt([i, j])
            # for t in TEMP_TARGETS:
            #     if t.inside(pt):
            #         target_data.set(pt_id, t.value)
            value = TEMP_TARGET.interp(pt)
            target_data.set(pt_id, value)
    target = cocoa.FieldMED()
    target.init("target", p, "on_nodes")
    target.set_values(target_data.data)
    target.init_io(p.output_prefix)
    target.print_vtk(0.0, 0)

    # discretized control
    control_data = cocoa.VectorFD(p.mesh.n_pts, 0.0)
    for j in range(p.mesh.n[1]):
        for i in range(p.mesh.n[0]):
            pt = p.mesh.pt([i, j])
            for c in CONTROL_REGIONS:
                if c.inside(pt):
                    control_data.set(i + p.mesh.n[0] * j, 1.0)
    control_mask = cocoa.FieldMED()
    control_mask.init("control_mask", p, "on_nodes")
    control_mask.set_values(control_data.data)
    control_mask.init_io(p.output_prefix)
    control_mask.print_vtk(0.0, 0)
    # p.set_field("control_mask", control_mask)

    # bcs
    p.bcs = [
        # forward
        cocoa.FDBCList2D(
            left=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.left,
                type=cocoa.FD_BC_TYPE.dirichlet,
                value=0.0,
                size=p.mesh.n[1],
            ),
            right=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.right,
                type=cocoa.FD_BC_TYPE.neumann,
                value=FLUX_EXT,
                size=p.mesh.n[1],
            ),
            bottom=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.bottom,
                type=cocoa.FD_BC_TYPE.neumann,
                value=FLUX_EXT,
                size=p.mesh.n[0],
            ),
            top=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.top,
                type=cocoa.FD_BC_TYPE.neumann,
                value=FLUX_EXT,
                size=p.mesh.n[0],
            ),
        ),
        # adjoint
        cocoa.FDBCList2D(
            left=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.left,
                type=cocoa.FD_BC_TYPE.dirichlet,
                value=0.0,
                size=p.mesh.n[1],
            ),
            right=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.right,
                type=cocoa.FD_BC_TYPE.neumann,
                value=0.0,
                size=p.mesh.n[1],
            ),
            bottom=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.bottom,
                type=cocoa.FD_BC_TYPE.neumann,
                value=0.0,
                size=p.mesh.n[0],
            ),
            top=cocoa.FDBC(
                side=cocoa.FD_BC_SIDE.top,
                type=cocoa.FD_BC_TYPE.neumann,
                value=0.0,
                size=p.mesh.n[0],
            ),
        ),
    ]

    # time
    p.time = 0.0
    p.final_time = FINAL_TIME
    p.dt = DT

    # la
    p.solver_type = cocoa.FD_SOLVER_TYPE.gauss_seidel
    p.m.init(size, 6)
    p.rhs.resize(size)

    # assembly
    def assembly(p: cocoa.ProblemFD2D):
        print("custom assembly heat oc")
        nx, ny = p.mesh.n

        for j in range(ny):
            for i in range(nx):
                id_f = i + j * p.mesh.n[0]
                id_a = i + j * p.mesh.n[0] + p.mesh.n_pts
                pt = p.mesh.pt([i, j])
                hx, hy = p.mesh.h

                # forward problem
                #   diagonal
                p.m.add(
                    id_f,
                    id_f,
                    1.0 + 2.0 * ALPHA * p.dt * (1.0 / (hx * hx) + 1.0 / (hy * hy)),
                )

                #   bottom
                value_bottom = -ALPHA * p.dt / hy**2
                if j > 0:
                    p.m.add(id_f, id_f - nx, value_bottom)
                else:
                    p.m.add(id_f, id_f + nx, value_bottom)
                    p.bcs[0].bottom.ghost_values.set(i, value_bottom)

                #   right
                value_right = -ALPHA * p.dt / hx**2
                if i < nx - 1:
                    p.m.add(id_f, id_f + 1, value_right)
                else:
                    p.m.add(id_f, id_f - 1, value_right)
                    p.bcs[0].right.ghost_values.set(j, value_right)

                #   top
                value_top = -ALPHA * p.dt / hy**2
                if j < ny - 1:
                    p.m.add(id_f, id_f + nx, value_top)
                else:
                    p.m.add(id_f, id_f - nx, value_top)
                    p.bcs[0].top.ghost_values.set(i, value_top)

                #   left
                value_left = -ALPHA * p.dt / hx**2
                if i > 0:
                    p.m.add(id_f, id_f - 1, value_left)
                else:
                    p.m.add(id_f, id_f + 1, value_left)
                    p.bcs[0].left.ghost_values.set(j, value_left)

                #   coupling with adjoint
                for c in CONTROL_REGIONS:
                    if c.inside(pt):
                        p.m.add(id_f, id_a, p.dt / BETA)

                #   rhs
                p.rhs.add(id_f, p.u_old[id_f])

                # adjoint problem
                #   diagonal
                p.m.add(
                    id_a,
                    id_a,
                    1.0 + 2.0 * ALPHA * p.dt * (1.0 / (hx * hx) + 1.0 / (hy * hy)),
                )

                #   bottom
                value_bottom = -ALPHA * p.dt / hy**2
                if j > 0:
                    p.m.add(id_a, id_a - nx, value_bottom)
                else:
                    p.m.add(id_a, id_a + nx, value_bottom)
                    p.bcs[1].bottom.ghost_values.set(i, value_bottom)

                #   right
                value_right = -ALPHA * p.dt / hx**2
                if i < nx - 1:
                    p.m.add(id_a, id_a + 1, value_right)
                else:
                    p.m.add(id_a, id_a - 1, value_right)
                    p.bcs[1].right.ghost_values.set(j, value_right)

                #   top
                value_top = -ALPHA * p.dt / hy**2
                if j < ny - 1:
                    p.m.add(id_a, id_a + nx, value_top)
                else:
                    p.m.add(id_a, id_a - nx, value_top)
                    p.bcs[1].top.ghost_values.set(i, value_top)

                #   left
                value_left = -ALPHA * p.dt / hx**2
                if i > 0:
                    p.m.add(id_a, id_a - 1, value_left)
                else:
                    p.m.add(id_a, id_a + 1, value_left)
                    p.bcs[1].left.ghost_values.set(j, value_left)

                #   coupling with forward problem
                # for t in TEMP_TARGETS:
                #     if t.inside(pt):
                #         p.m.add(id_a, id_f, -p.dt)
                #         p.rhs.add(id_a, -t.value * p.dt)
                if target.at(id_f) > 1e-12:
                    p.m.add(id_a, id_f, -p.dt)
                    p.rhs.add(id_a, -target.at(id_f) * p.dt)

                #   rhs
                p.rhs.add(id_a, p.u_old[id_a])

        p.m.close()

    p.set_custom_assembly(assembly)

    # p.eqn_type = cocoa.EQN_TYPE.heat_oc

    print("assemblies: {}", p.assemblies)


def setup(p: cocoa.ProblemFD2D):
    p.name = "fd2d_custom"
    p.debug = False

    # mesh
    p.mesh = cocoa.MeshFD2D((0.0, 0.0), (1.0, 1.0), (10, 10))

    # coupling
    p.coupling_type = cocoa.COUPLING_TYPE.medcoupling

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


def compute_functional_distributed(p: cocoa.ProblemFD2D) -> float:
    func = 0.0

    da = p.mesh.h[0] * p.mesh.h[1]

    for j in range(p.mesh.n[1]):
        for i in range(p.mesh.n[0]):
            id_f = i + p.mesh.n[0] * j
            id_a = id_f + p.mesh.n_pts
            pt = p.mesh.pt([i, j])

            # for t in TEMP_TARGETS:
            #     if t.inside(pt):
            #         func += 0.5 * (p.u[id_f] - t.value) * (p.u[id_f] - t.value) * da
            value = TEMP_TARGET.interp(pt)
            func += 0.5 * (p.u[id_f] - value) * (p.u[id_f] - value) * da

        for c in CONTROL_REGIONS:
            if c.inside(pt):
                func += 0.5 * p.u[id_a] * p.u[id_a] * da / BETA

    return func


def main():
    p = cocoa.ProblemFD2D()

    # p.setup(config_file="fd2d_heat.dat")
    # setup(p)

    # p.setup(config_file="fd2d_oc.dat")
    setup_oc(p)

    p.print()

    with open(p.output_prefix / "func.csv", "w") as f:
        f.write("time,num_iters,func\n")
        func_old = 1.0
        prepareToBreak = False
        while p.run():
            p.advance()
            num_iters = p.solve()
            p.print()

            func = compute_functional_distributed(p)
            print(f"functional: {func:.12e}")
            f.write(f"{p.time:.8e},{num_iters:8d},{func:.16e}\n")

            if np.fabs(func - func_old) < FUNC_TOL:
                if prepareToBreak:
                    break
                else:
                    prepareToBreak = True

            func_old = func

            if num_iters == 0:
                if prepareToBreak:
                    break
                else:
                    prepareToBreak = True


if __name__ == "__main__":
    main()
