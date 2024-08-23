import pycocoa

pycocoa.setFD1DAssemblies()

p1 = pycocoa.ProblemFD1D()

p1.setup({"config_file": "fd1d1.dat"})
p1.print()

while p1.run():
    p1.advance()
    p1.solve()
    p1.print()

