import cocoa

p1 = cocoa.ProblemFD1D()
p2 = cocoa.ProblemFD1D()

p1.setup(config_file="fd1d_heat.dat")
p2.setup(config_file="fd1d_hc.dat")

p1.print()
p2.print()

# c = cocoa.CouplingSimple(cocoa.COUPLING_SCOPE.volume)
c = cocoa.CouplingSimple(cocoa.COUPLING_SCOPE.volume)
interface_src = cocoa.CouplingInterface(p1, ["T"])
interface_tgt = cocoa.CouplingInterface(p2, ["Tcfd"])
c.setup(interface_src, interface_tgt, cocoa.INTERPOLATION_METHOD.p1p1)

while p1.run() or p2.run():
    p1.advance()
    p1.solve()
    p1.print()

    c.project(name_src="T", name_tgt="Tcfd")
    p2.advance()
    p2.solve()
    p2.print()

    c.print_vtk()
