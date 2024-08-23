import pycocoa

p1 = pycocoa.ProblemProXPDE()
p1.setup(config_file="proxpde1.yaml")
p1.print()

p2 = pycocoa.ProblemProXPDE()
p2.setup(config_file="proxpde2.yaml")
p2.print()

print(f"coupling type: {p1.coupling_type}")
assert p1.coupling_type == p2.coupling_type

c = pycocoa.CouplingMED()
c.setup(problem_src=p1, problem_tgt=p2)

while p1.run() or p2.run():
    p1.advance()
    p1.solve()
    p1.print()

    c.project(name_src="T", name_tgt="Tsys")
    p2.advance()
    p2.solve()
    p2.print()
