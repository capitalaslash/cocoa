import pycocoa

p = pycocoa.ProblemProXPDEHeat()

p.setup(config_file="proxpde1.yaml")
p.print()

while p.run():
    p.advance()
    p.solve()
    p.print()
