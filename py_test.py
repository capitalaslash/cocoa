import fd1d

p = fd1d.ProblemFD1D()

p.setup({"config_file": "fd1d.dat"})

while(p.run()):
    p.advance()
    p.solve()
    p.print()

