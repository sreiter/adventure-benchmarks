-- Copyright (c) 2010-2017:  G-CSC, Goethe University Frankfurt
-- Authors: Andreas Vogel, Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")
ug_load_script("util/profiler_util.lua")

-- Parse parameters and print help
dim			= 3
gridName	= util.GetParam("-grid", "box-tet.ugx", "filename of underlying grid")
numRefs		= util.GetParamNumber("-numRefs", 3, "number of refinements")
writeSolution	= util.HasParamOption("-writeSolution", "If specified, the solution will be written to a file.")
redistElemThreshold = util.GetParamNumber("-redistElemThreshold", 8, "Number of elements that each target process should receive")
redistProcs = util.GetParamNumber("-redistProcs", 64, "Number of processes to which each process redistributes on a distribution level.")

util.CheckAndPrintHelp("Laplace-HPC");
util.PrintArguments()


-- initialize ug with the world dimension and the algebra type
InitUG(dim, AlgebraType("CPU", 1));


-- Load a domain without initial refinements.
requiredSubsets = {"inner", "bndTop", "bndBottom"}
dom = util.CreateDomain(gridName, 0, requiredSubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")

-- This balancing setup makes sense for structured grids with uniform refinement
balancerDesc = {
		partitioner = {
			name = "staticBisection",
			clusteredSiblings = false
		},

		hierarchy = {
			name 						= "noRedists",
			minElemsPerProcPerLevel		= redistElemThreshold,
			maxRedistProcs				= redistProcs,
		},
	}

util.refinement.CreateRegularHierarchy(dom, numRefs, true, balancerDesc)


-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


-- set up discretization
elemDisc = ConvectionDiffusion("c", "inner", "fv1")
elemDisc:set_diffusion(1.0)

dirichletBND = DirichletBoundary()
dirichletBND:add(1, "c", "bndTop")
dirichletBND:add(0, "c", "bndBottom")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)


-- set up solver (using 'util/solver_util.lua')
solverDesc = {
	type = "bicgstab",
	precond = {
		type		= "gmg",
		approxSpace	= approxSpace,
		smoother	= "jac",
		baseSolver	= "lu"
	},

--	Settings make sure that exactly 10 iterations are performed.
--	This is important to compare solver run times.
--	NOTE: This means that the solution step will not be considered successful
	convCheck = {
		name		= "standard",
		iterations	= 10,
		absolute	= 1e-24,
		reduction	= 1e-24,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)


print("\nsolving...")
A = AssembledLinearOperator(domainDisc)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)


solver:init(A, u)
solver:apply(u, b)

print("NOTE: The hard limit on iteration numbers is intentional here, to allow")
print("      for a better comparability of run-times. Please make sure that the")
print("      reached defect only varies slightly regarding different refinement")
print("      and process numbers.")

if writeSolution == true then
	solFileName = "sol_laplace_"..dim.."d"
	print("writing solution to '" .. solFileName .. "'...")
	WriteGridFunctionToVTK(u, solFileName)
end

if GetProfilerAvailable() == true then
	WriteProfileData("laplace.pdxml")
	print("")
	print("--- PROFILES --- ")
	util.PrintProfile_TotalTime("main                           ")
	util.PrintProfile_TotalTime("  assemble_linear              ")
	util.PrintProfile_TotalTime("  LS_InitPrecond               ")
	util.PrintProfile_TotalTime("  LS_ApplyReturnDefect         ")
end

print("done")
