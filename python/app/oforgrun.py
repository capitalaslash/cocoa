import cocoa

p = cocoa.ProblemOForg()

p.setup(case_dir="oforg1")
p.print()

while p.run():
    p.advance()
    p.solve()
    p.print()
