import fd1d

fd1d.setFD1DAssemblies()

p1 = fd1d.ProblemFD1D()

p1.setup({"config_file": "fd1d1.dat"})
p1.print()

while p1.run():
    p1.advance()
    p1.solve()
    p1.print()
