import cocoa  # type: ignore

p = cocoa.ProblemFD1D()

p.setup(config_file="fd1d_heat.dat")
p.print()

while p.run():
    p.advance()
    p.solve()
    p.print()
