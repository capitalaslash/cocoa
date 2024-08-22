import fd1d
import coupling_simple as cs

fd1d.setFD1DAssemblies()

p1 = fd1d.ProblemFD1D()
p1.setup({"config_file": "fd1d1.dat"})
p1.print()

p2 = fd1d.ProblemFD1D()
p2.setup({"config_file": "fd1d2.dat"})
p2.print()

c = cs.CouplingSimple()
c.setup(p1, p2)

while p1.run() or p2.run():
    p1.advance()
    p1.solve()
    p1.print()

    c.project("u", "uExternal")
    p2.advance()
    p2.solve()
    p2.print()
