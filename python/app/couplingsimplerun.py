import pycocoa

p1 = pycocoa.ProblemFD1D()
p1.setup(config_file="fd1d_heat.dat")
p1.print()

p2 = pycocoa.ProblemFD1D()
p2.setup(config_file="fd1d_hc.dat")
p2.print()

c = pycocoa.CouplingSimple()
# c = pycocoa.CouplingMED()
c.setup(problem_src=p1, problem_tgt=p2)

while p1.run() or p2.run():
    p1.advance()
    p1.solve()
    p1.print()

    c.project(name_src="u", name_tgt="uExternal")
    p2.advance()
    p2.solve()
    p2.print()
