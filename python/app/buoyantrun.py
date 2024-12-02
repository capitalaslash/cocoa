import pycocoa

pNS = pycocoa.ProblemProXPDENS()
pNS.setup(config_file="proxpde_buoyant_ns.yaml")

pHeat = pycocoa.ProblemProXPDEHeat()
pHeat.setup(config_file="proxpde_buoyant_heat.yaml")

couplingHeatToNS = pycocoa.CouplingMED()
couplingHeatToNS.setup(problem_src=pHeat, problem_tgt=pNS)

couplingNSToHeat = pycocoa.CouplingMED()
couplingNSToHeat.setup(problem_src=pNS, problem_tgt=pHeat)

pNS.print()
pHeat.print()

while pNS.run() or pHeat.run():
    couplingHeatToNS.project("T")
    pNS.advance()
    pNS.solve()
    pNS.print()

    couplingNSToHeat.project("vel")
    pHeat.advance()
    pHeat.solve()
    pHeat.print()
