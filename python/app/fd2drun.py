import pycocoa

p = pycocoa.ProblemFD2D()

p.setup(config_file="fd2d_heat.dat")
p.print()

while p.run():
    p.advance()
    p.solve()
    p.print()
