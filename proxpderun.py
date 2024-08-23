import pycocoa

p1 = pycocoa.ProblemProXPDE()

p1.setup(config_file="proxpde1.yaml")
p1.print()

while p1.run():
    p1.advance()
    p1.solve()
    p1.print()
