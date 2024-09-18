import pycocoa

p = pycocoa.ProblemProXPDEHeat()

p.setup(config_file="proxpde_heat.yaml")
p.print()

while p.run():
    p.advance()
    p.solve()
    p.print()
