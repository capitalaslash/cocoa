import cocoa

pNS = cocoa.ProblemProXPDENS()
pNS.setup(config_file="proxpde_buoyant_ns.yaml")

pHeat = cocoa.ProblemProXPDEHeat()
pHeat.setup(config_file="proxpde_buoyant_heat.yaml")

couplingHeatToNS = cocoa.CouplingMED()
couplingHeatToNS.setup(problem_src=pHeat, problem_tgt=pNS)

couplingNSToHeat = cocoa.CouplingMED()
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
