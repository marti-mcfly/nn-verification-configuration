import os
import subprocess
import argparse
import sys
import shutil
import random
import socket
import string	

def check_arg(args=None):
	parser = argparse.ArgumentParser(description='Script to run pipeline')
	parser.add_argument('--name', help='Name of MIP instance file', required='True')

	# parameters for Gurobi..
	parser.add_argument('--auto_threads_on', default=0, type=int)
	parser.add_argument('--threads', default=1, type=int)
	
	parser.add_argument('--auto_issmethod_on', default=1, type=int)
	parser.add_argument('--issmethod', default=-1, type=int)
	
	parser.add_argument('--auto_presolve_on', default=1, type=int)
	parser.add_argument('--presolve', default=-1, type=int)
	
	parser.add_argument('--auto_aggfill_on', default=1, type=int)
	parser.add_argument('--aggfill', default=-1, type=int)
	
	parser.add_argument('--aggregate', default=1, type=int)
	parser.add_argument('--dualreductions', default=1, type=int)
	parser.add_argument('--precrush', default=0, type=int)
	parser.add_argument('--predeprow', default=1, type=int)
	
	parser.add_argument('--auto_predual_on', default=1, type=int)
	parser.add_argument('--predual', default=-1, type=int)
	
	parser.add_argument('--auto_prepasses_on', default=1, type=int)
	parser.add_argument('--prepasses', default=-1, type=int)
	
	parser.add_argument('--auto_presparsify_on', default=1, type=int)
	parser.add_argument('--presparsify', default=-1, type=int)
	
	parser.add_argument('--crossoverbasis', default=0, type=int)
	
	parser.add_argument('--auto_barcorrectors', default=1, type=int)
	parser.add_argument('--barcorrectors', default=-1, type=int)
	
	parser.add_argument('--auto_barhomogeneous', default=1, type=int)
	parser.add_argument('--barhomogeneous', default=-1, type=int)
	
	parser.add_argument('--auto_barorder', default=1, type=int)
	parser.add_argument('--barorder', default=-1, type=int)
	
	parser.add_argument('--auto_crossover', default=1, type=int)
	parser.add_argument('--crossover', default=-1, type=int)
	
	parser.add_argument('--infunbdinfo', default=0, type=int)
	parser.add_argument('--perturbvalue', default=0.0002, type=float)
	
	parser.add_argument('--auto_normadjust_on', default=1, type=int)
	parser.add_argument('--normadjust', default=-1, type=int)
	
	parser.add_argument('--auto_quad_on', default=1, type=int)
	parser.add_argument('--quad', default=-1, type=int)
	
	parser.add_argument('--auto_sifting_on', default=1, type=int)
	parser.add_argument('--sifting', default=-1, type=int)
	
	parser.add_argument('--auto_siftmethod_on', default=1, type=int)
	parser.add_argument('--siftmethod', default=-1, type=int)
	
	parser.add_argument('--auto_simplexpricing_on', default=1, type=int)
	parser.add_argument('--simplexpricing', default=-1, type=int)
	
	parser.add_argument('--heuristics', default=0.05, type=float)
	parser.add_argument('--improvestartgap', default=0.0, type=float)
	parser.add_argument('--improvestartnodes', default=0.0, type=float)
	parser.add_argument('--improvestarttime', default=0.0, type=float)
	parser.add_argument('--mipfocus', default=0, type=int)
	parser.add_argument('--submipnodes', default=500, type=int)
	parser.add_argument('--partitionplace', default=15, type=int)
	
	parser.add_argument('--auto_branchdir_on', default=1, type=int)
	parser.add_argument('--branchdir', default=0, type=int)
	
	parser.add_argument('--auto_degenmoves_on', default=1, type=int)
	parser.add_argument('--degenmoves', default=-1, type=int)
	
	parser.add_argument('--auto_disconnected_on', default=1, type=int)
	parser.add_argument('--disconnected', default=-1, type=int)
	
	parser.add_argument('--auto_minrelnodes_on', default=1, type=int)
	parser.add_argument('--minrelnodes', default=-1, type=int)
	
	parser.add_argument('--auto_nodemethod_on', default=1, type=int)
	parser.add_argument('--nodemethod', default=-1, type=int)
	
	parser.add_argument('--auto_pumppasses_on', default=1, type=int)
	parser.add_argument('--pumppasses', default=-1, type=int)
	
	parser.add_argument('--auto_rins_on', default=1, type=int)
	parser.add_argument('--rins', default=-1, type=int)
	
	parser.add_argument('--shut_off_mip_start_processing', default=0, type=int) 
	parser.add_argument('--auto_startnodelimit', default=1, type=int)
	parser.add_argument('--startnodelimit', default=-2, type=int) #default can be -1 (will be set later)
	
	parser.add_argument('--auto_symmetry_on', default=1, type=int)
	parser.add_argument('--symmetry', default=-1, type=int)
	
	parser.add_argument('--auto_varbranch_on', default=1, type=int)
	parser.add_argument('--varbranch', default=-1, type=int)
	
	parser.add_argument('--auto_zeroobjnodes_on', default=1, type=int)
	parser.add_argument('--zeroobjnodes', default=-1, type=int)
	
	parser.add_argument('--auto_cuts_on', default=1, type=int)
	parser.add_argument('--cuts', default=-1, type=int)
	
	parser.add_argument('--auto_bqpcuts_on', default=1, type=int)
	parser.add_argument('--bqpcuts', default=-1, type=int)
	
	parser.add_argument('--auto_cliquecuts_on', default=1, type=int)
	parser.add_argument('--cliquecuts', default=-1, type=int)
	
	parser.add_argument('--auto_covercuts_on', default=1, type=int)
	parser.add_argument('--covercuts', default=-1, type=int)
	
	parser.add_argument('--auto_cutaggpasses_on', default=1, type=int)
	parser.add_argument('--cutaggpasses', default=-1, type=int)
	
	parser.add_argument('--auto_cutpasses_on', default=1, type=int)
	parser.add_argument('--cutpasses', default=-1, type=int)
	
	parser.add_argument('--auto_flowcovercuts_on', default=1, type=int)
	parser.add_argument('--flowcovercuts', default=-1, type=int)
	
	parser.add_argument('--auto_flowpathcuts_on', default=1, type=int)
	parser.add_argument('--flowpathcuts', default=-1, type=int)
	
	parser.add_argument('--auto_gomorypasses_on', default=1, type=int)
	parser.add_argument('--gomorypasses', default=-1, type=int)
	
	parser.add_argument('--auto_gubcovercuts_on', default=1, type=int)
	parser.add_argument('--gubcovercuts', default=-1, type=int)
	
	parser.add_argument('--auto_impliedcuts_on', default=1, type=int)
	parser.add_argument('--impliedcuts', default=-1, type=int)
	
	parser.add_argument('--auto_infproofcuts_on', default=1, type=int)
	parser.add_argument('--infproofcuts', default=-1, type=int)
	
	parser.add_argument('--auto_mipsepcuts_on', default=1, type=int)
	parser.add_argument('--mipsepcuts', default=-1, type=int)
	
	parser.add_argument('--auto_mircuts_on', default=1, type=int)
	parser.add_argument('--mircuts', default=-1, type=int)
	
	parser.add_argument('--auto_modkcuts_on', default=1, type=int)
	parser.add_argument('--modkcuts', default=-1, type=int)
	
	parser.add_argument('--auto_networkcuts_on', default=1, type=int)
	parser.add_argument('--networkcuts', default=-1, type=int)
	
	parser.add_argument('--auto_projimpliedcuts_on', default=1, type=int)
	parser.add_argument('--projimpliedcuts', default=-1, type=int)
	
	parser.add_argument('--auto_relaxliftcuts_on', default=1, type=int)
	parser.add_argument('--relaxliftcuts', default=-1, type=int)
	
	parser.add_argument('--auto_rltcuts_on', default=1, type=int)
	parser.add_argument('--rltcuts', default=-1, type=int)
	
	parser.add_argument('--auto_strongcgcuts_on', default=1, type=int)
	parser.add_argument('--strongcgcuts', default=-1, type=int)
	
	parser.add_argument('--auto_submipcuts_on', default=1, type=int)
	parser.add_argument('--submipcuts', default=-1, type=int)
	
	parser.add_argument('--auto_zerohalfcuts_on', default=1, type=int)
	parser.add_argument('--zerohalfcuts', default=-1, type=int)
	
	return parser.parse_args(args)

def set_parameters(model, args):
	#setting Gurobi parameters..........
	model.Params.Threads = args.threads
	model.Params.IISMethod = args.issmethod
	model.Params.Presolve = args.presolve
	model.Params.AggFill = args.aggfill
	model.Params.Aggregate = args.aggregate
	model.Params.DualReductions = args.dualreductions
	model.Params.PreCrush = args.precrush
	model.Params.PreDeProw = args.predeprow
	model.Params.PreDual = args.predual
	model.Params.PrePasses = args.prepasses
	model.Params.PreSparsify = args.presparsify
	model.Params.CrossoverBasis = args.crossoverbasis
	model.Params.BarCorrectors = args.barcorrectors
	model.Params.BarHomogeneous = args.barhomogeneous
	model.Params.BarOrder = args.barorder
	model.Params.Crossover = args.crossover
	model.Params.InfUnbdInfo = args.infunbdinfo
	model.Params.PerturbValue = args.perturbvalue
	model.Params.NormAdjust = args.normadjust
	model.Params.Quad = args.quad
	model.Params.Sifting = args.sifting
	model.Params.SiftMethod = args.siftmethod
	model.Params.SimplexPricing = args.simplexpricing
	model.Params.Heuristics = args.heuristics
	model.Params.ImproveStartGap = args.improvestartgap
	model.Params.ImproveStartNodes = args.improvestartnodes
	model.Params.ImproveStartTime = args.improvestarttime
	model.Params.MIPFocus = args.mipfocus
	model.Params.SubMIPNodes = args.submipnodes
	model.Params.PartitionPlace = args.partitionplace
	model.Params.BranchDir = args.branchdir
	model.Params.DegenMoves = args.degenmoves
	model.Params.Disconnected = args.disconnected
	model.Params.MinRelNodes = args.minrelnodes
	model.Params.NodeMethod = args.nodemethod
	model.Params.PumpPasses = args.pumppasses
	model.Params.RINS = args.rins
	if args.shut_off_mip_start_processing == 1 and args.auto_startnodelimit == 0:
		startnodelimit = -2
		model.Params.StartNodeLimit = startnodelimit
	elif args.shut_off_mip_start_processing == 0 and args.auto_startnodelimit == 1:
		startnodelimit = -1
		model.Params.StartNodeLimit = startnodelimit
	else:
		model.Params.StartNodeLimit = args.startnodelimit
	model.Params.Symmetry = args.symmetry
	model.Params.VarBranch = args.varbranch
	model.Params.ZeroObjNodes = args.zeroobjnodes
	model.Params.Cuts = args.cuts
	model.Params.BQPCuts = args.bqpcuts
	model.Params.CliqueCuts = args.cliquecuts
	model.Params.CoverCuts = args.covercuts
	model.Params.CutAggPasses = args.cutaggpasses
	model.Params.FlowCoverCuts = args.flowcovercuts
	model.Params.FlowPathCuts = args.flowpathcuts
	model.Params.GomoryPasses = args.gomorypasses
	model.Params.GubCoverCuts = args.gubcovercuts
	model.Params.ImpliedCuts = args.impliedcuts
	model.Params.INFProofCuts = args.infproofcuts
	model.Params.MIPSepCuts = args.mipsepcuts
	model.Params.MIRCuts = args.mircuts
	model.Params.ModKCuts = args.modkcuts
	model.Params.NetworkCuts = args.networkcuts
	model.Params.ProjImpliedCuts = args.projimpliedcuts
	model.Params.RelaxLiftCuts = args.relaxliftcuts
	model.Params.RLTCuts = args.rltcuts
	model.Params.StrongCGCuts = args.strongcgcuts
	model.Params.SubMIPCuts = args.submipcuts
	model.Params.ZeroHalfCuts = args.zerohalfcuts

def variable_order(reordering, symmetric, converging, window):
	if(symmetric == 1 and converging == 1):
		print('Variable order algorithm: CUDD_REORDER_SYMM_SIFT_CONV') 
		return 7
	elif(symmetric == 1):
		print('Variable order algorithm: CUDD_REORDER_SYMM_SIFT')
		return 6
	elif (converging == 1):
		if(reordering == 'sift'):
			print('Variable order algorithm: CUDD_REORDER_SIFT_CONVERGE')
			return 5
		elif(reordering == 'group-sift'):
			print('Variable order algorithm: CUDD_REORDER_GROUP_SIFT_CONV')
			return 15
		elif(reordering == 'window'):
			if(window == 2):
				print('Variable order algorithm: CUDD_REORDER_WINDOW2_CONV')
				return 11
			elif(window == 3):
				print('Variable order algorithm: CUDD_REORDER_WINDOW3_CONV')
				return 12
			elif(window == 4):
				print('Variable order algorithm: CUDD_REORDER_WINDOW4_CONV')	
				return 13
	else:
		if(reordering == 'sift'):
			print('Variable order algorithm: CUDD_REORDER_SIFT')
			return 4
		elif(reordering == 'group-sift'):
			print('Variable order algorithm: CUDD_REORDER_GROUP_SIFT')
			return 14
		elif(reordering == 'window'):
			if(window == 2):
				print('Variable order algorithm: CUDD_REORDER_WINDOW2')
				return 8
			elif(window == 3):
				print('Variable order algorithm: CUDD_REORDER_WINDOW3')
				return 9
			elif(window == 4):
				print('Variable order algorithm: CUDD_REORDER_WINDOW4')	
				return 10
		elif(reordering == 'genetic'):
			print('Variable order algorithm: CUDD_REORDER_GENETIC')
			return 17
		elif(reordering == 'annealing'):
			print('Variable order algorithm: CUDD_REORDER_ANNEALING')
			return 18
		elif(reordering == 'random'):
			print('Variable order algorithm: CUDD_REORDER_RANDOM')
			return 2
			
args = check_arg(sys.argv[1:])
reorder_type = 14 #default
dyn_reorder_type = 14

print('Name of MIP instance file =', args.name)

from gurobipy import *

model = read(args.name)

set_parameters(model, args)

print('Solving..')
sys.stdout.flush()

model.optimize()

# print("RUN SUCCESSFULLY COMPLETED")
