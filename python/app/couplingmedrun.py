import math

import cocoa


def source(p):
    return math.sin(math.pi * p[1])


p1 = cocoa.ProblemProXPDEHeat()
p1.setup(config_file="proxpde_heat.yaml")
p1.set_source(source)

p2 = cocoa.ProblemProXPDEHeat()
p2.setup(config_file="proxpde_hc.yaml")

c = cocoa.CouplingMED()
c.setup(problem_src=p1, problem_tgt=p2)

p1.print()
p2.print()

while p1.run() or p2.run():
    p1.advance()
    p1.solve()
    p1.print()

    c.project(name="Tcfd")
    p2.advance()
    p2.solve()
    p2.print()
