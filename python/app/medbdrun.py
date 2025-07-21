import cocoa

p1 = cocoa.ProblemProXPDENS()
p2 = cocoa.ProblemProXPDENS()

p1.setup(config_file="proxpde_med_bd1.yaml")
p2.setup(config_file="proxpde_med_bd2.yaml")

p1.print()
p2.print()

coupling12 = cocoa.CouplingMED(cocoa.COUPLING_SCOPE.boundary)
coupling12.setup(
    interface_src=cocoa.CouplingInterface(p1, "top", ["vel"]),
    interface_tgt=cocoa.CouplingInterface(p2, "bottom", ["vel"]),
    method=cocoa.INTERPOLATION_METHOD.p1p1,
)
coupling21 = cocoa.CouplingMED(cocoa.COUPLING_SCOPE.boundary)
coupling21.setup(
    interface_src=cocoa.CouplingInterface(p2, "bottom", ["p"]),
    interface_tgt=cocoa.CouplingInterface(p1, "top", ["p"]),
    method=cocoa.INTERPOLATION_METHOD.p1p1,
)

while p1.run() | p2.run():
    p1.advance()
    p2.advance()

    for _ in range(2):
        coupling21.project("p")
        p1.solve()
        coupling12.project("vel")
        p2.solve()

    p1.print()
    p2.print()
