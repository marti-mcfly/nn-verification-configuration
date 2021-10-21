import sys
sys.path.append('.')

import argparse
import numpy as np
import pickle
import random
from timeit import default_timer as timer

#### to stop seeing various warning and system messages
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
import keras
sys.stderr = stderr

import logging
import tensorflow as tf
tf.get_logger().setLevel(logging.ERROR)
####

# from resources.acas.acasprop import acas_properties, acas_denormalise_input, acas_normalise_input
from src.parameters import SplittingParam, Param
from src.specifications import GenSpec, LRob
from src.venusverifier import VenusVerifier

verifier_result_to_string = {"True": "NOT Satis", "False": "Satisfied", "Interrupted": "Interrupt", "Timeout": "Timed Out"}

def verify_acas_property(options):
    random.seed(options.acas_prop)

    prop = acas_properties[options.acas_prop]

    input_bounds = (np.array(prop['bounds']['lower']), np.array(prop['bounds']['upper']))
    spec = GenSpec(input_bounds, prop['output'])

    encoder_params = Param()
    encoder_params.TIME_LIMIT = options.timeout
    encoder_params.XLAYER_DEP_CONSTRS = options.offline_dep
    encoder_params.LAYER_DEP_CONSTRS = options.offline_dep
    encoder_params.XLAYER_DEP_CUTS = options.online_dep
    encoder_params.LAYER_DEP_CUTS = options.online_dep
    encoder_params.IDEAL_CUTS = options.ideal_cuts
    encoder_params.WORKERS_NUMBER = options.workers

    splitting_params = SplittingParam()
    splitting_params.FIXED_RATIO_CUTOFF = options.st_ratio
    splitting_params.DEPTH_POWER = options.depth_power
    splitting_params.STARTING_DEPTH_PARALLEL_SPLITTING_PROCESSES = options.splitters
    
    verifier = VenusVerifier(options.net, spec, encoder_params, splitting_params, options.print)

    start = timer()
    result, job_id, extra = verifier.verify()
    end = timer()
    runtime = end - start

    result = verifier_result_to_string[result]

    print("{} over {}".format(prop['name'], options.net),
          "is", result, "in {:9.4f}s".format(runtime), "job n", job_id)
    if result == True:
        denormalised_ctx = acas_denormalise_input(extra)
        ctx = extra.reshape(1, -1)
        network_output = nmodel.predict(x=ctx, batch_size=1)
        print("\t\tCounter-example:", list(extra))
        print("\t\tDenormalised   :", list(denormalised_ctx))
        print("\t\tNetwork output :", list(network_output[0]))
    print("")

    return result, runtime, job_id


def verify_local_robustness(options):
    # load image
    with open(options.lrob_input, 'rb') as pickle_file:
        data = pickle.load(pickle_file)
    image = data[0]
    label = data[1]
    num_classes = 10
    # range of pixel values
    min_value = 0
    max_value = 1

    # create specification
    spec = LRob(image, label, num_classes, options.lrob_radius, min_value, max_value)

    # set verifier's parameters
    encoder_params = Param()
    encoder_params.TIME_LIMIT = options.timeout
    encoder_params.THREADS = options.threads
    encoder_params.AUTO_ISSMETHOD_ON = options.auto_issmethod_on
    encoder_params.ISSMETHOD = options.issmethod
    encoder_params.AUTO_PRESOLVE_ON = options.auto_presolve_on 
    encoder_params.PRESOLVE = options.presolve
    encoder_params.AUTO_AGGFILL_ON = options.auto_aggfill_on
    encoder_params.AGGFILL = options.aggfill
    encoder_params.AGGREGATE = options.aggregate 
    encoder_params.DUALREDUCTIONS = options.dualreductions 
    encoder_params.PRECRUSH = options.precrush
    encoder_params.PREDEPROW = options.predeprow 
    encoder_params.AUTO_PREDUAL_ON = options.auto_predual_on 
    encoder_params.PREDUAL = options.predual 
    encoder_params.AUTO_PREPASSES_ON = options.auto_prepasses_on
    encoder_params.PREPASSES = options.prepasses
    encoder_params.AUTO_PRESPARSIFY_ON = options.auto_presparsify_on 
    encoder_params.PRESPARSIFY = options.presparsify
    encoder_params.CROSSOVERBASIS = options.crossoverbasis 
    encoder_params.AUTO_BARCORRECTORS = options.auto_barcorrectors 
    encoder_params.BARCORRECTORS = options.barcorrectors 
    encoder_params.AUTO_BARHOMOGENEOUS = options.auto_barhomogeneous 
    encoder_params.BARHOMOGENEOUS = options.barhomogeneous 
    encoder_params.AUTO_BARORDER = options.auto_barorder 
    encoder_params.BARORDER = options.barorder 
    encoder_params.AUTO_CROSSOVER = options.auto_crossover 
    encoder_params.CROSSOVER = options.crossover
    encoder_params.INFUNBDINFO = options.infunbdinfo 
    encoder_params.PERTURBVALUE = options.perturbvalue 
    encoder_params.AUTO_NORMADJUST_ON = options.normadjust 
    encoder_params.AUTO_QUAD_ON = options.auto_quad_on 
    encoder_params.QUAD = options.quad 
    encoder_params.AUTO_SIFTING_ON = options.auto_sifting_on 
    encoder_params.SIFTING = options.sifting 
    encoder_params.AUTO_SIFTMETHOD_ON = options.auto_siftmethod_on 
    encoder_params.SIFTMETHOD = options.siftmethod 
    encoder_params.AUTO_SIMPLEXPRICING_ON = options.auto_simplexpricing_on 
    encoder_params.SIMPLEXPRICING = options.simplexpricing 
    encoder_params.HEURISTICS = options.heuristics 
    encoder_params.IMPROVESTARTGAP = options.improvestartgap 
    encoder_params.IMPROVESTARTNODES = options.improvestartnodes 
    encoder_params.IMPROVESTARTTIME = options.improvestarttime 
    encoder_params.MIPFOCUS = options.mipfocus 
    encoder_params.SUBMIPNODES = options.submipnodes 
    encoder_params.PARTITIONPLACE = options.partitionplace 
    encoder_params.AUTO_BRANCHDIR_ON = options.auto_branchdir_on
    encoder_params.BRANCHDIR = options.branchdir 
    encoder_params.AUTO_DEGENMOVES_ON = options.auto_degenmoves_on 
    encoder_params.DEGENMOVES = options.degenmoves 
    encoder_params.AUTO_DISCONNECTED_ON = options.auto_disconnected_on 
    encoder_params.DISCONNECTED = options.disconnected
    encoder_params.AUTO_MINRELNODES_ON = options.auto_minrelnodes_on 
    encoder_params.MINRELNODES = options.minrelnodes 
    encoder_params.AUTO_NODEMETHOD_ON = options.auto_nodemethod_on
    encoder_params.NODEMETHOD = options.nodemethod 
    encoder_params.AUTO_PUMPPASSES_ON = options.pumppasses 
    encoder_params.PUMPPASSES = options.pumppasses 
    encoder_params.AUTO_RINS_ON = options.auto_rins_on 
    encoder_params.RINS = options.rins
    encoder_params.AUTO_STARTNODELIMIT = options.auto_startnodelimit 
    encoder_params.STARTNODELIMIT = options.startnodelimit 
    encoder_params.AUTO_SYMMETRY_ON = options.auto_symmetry_on 
    encoder_params.SYMMETRY = options.symmetry 
    encoder_params.AUTO_VARBRANCH_ON = options.auto_varbranch_on 
    encoder_params.VARBRANCH = options.varbranch 
    encoder_params.AUTO_ZEROOBJNODES_ON = options.auto_zeroobjnodes_on 
    encoder_params.ZEROOBJNODES = options.zeroobjnodes  
    encoder_params.AUTO_CUTS_ON = options.auto_cuts_on 
    encoder_params.CUTS = options.cuts 
    encoder_params.AUTO_BQPCUTS_ON = options.auto_bqpcuts_on
    encoder_params.BQPCUTS = options.bqpcuts 
    encoder_params.AUTO_CLIQUECUTS_ON = options.auto_cliquecuts_on 
    encoder_params.CLIQUECUTS = options.cliquecuts 
    encoder_params.AUTO_COVERCUTS_ON = options.auto_covercuts_on 
    encoder_params.COVERCUTS = options.covercuts 
    encoder_params.AUTO_CUTAGGPASSES_ON = options.auto_cutaggpasses_on 
    encoder_params.CUTAGGPASSES = options.cutaggpasses 
    encoder_params.AUTO_FLOWCOVERCUTS_ON = options.auto_flowcovercuts_on 
    encoder_params.FLOWCOVERCUTS = options.flowcovercuts 
    encoder_params.AUTO_FLOWPATHCUTS_ON = options.auto_flowpathcuts_on 
    encoder_params.FLOWPATHCUTS = options.flowpathcuts 
    encoder_params.AUTO_GOMORYPASSES_ON = options.auto_gomorypasses_on 
    encoder_params.GOMORYPASSES = options.gomorypasses 
    encoder_params.AUTO_GUBCOVERCUTS_ON = options.auto_gubcovercuts_on 
    encoder_params.GUBCOVERCUTS = options.gubcovercuts 
    encoder_params.AUTO_IMPLIEDCUTS_ON = options.auto_impliedcuts_on 
    encoder_params.IMPLIEDCUTS = options.impliedcuts 
    encoder_params.AUTO_INFPROOFCUTS_ON = options.auto_infproofcuts_on 
    encoder_params.INFPROOFCUTS = options.infproofcuts 
    encoder_params.AUTO_MIPSEPCUTS_ON = options.auto_mipsepcuts_on 
    encoder_params.MIPSEPCUTS = options.mipsepcuts 
    encoder_params.AUTO_MIRCUTS_ON = options.auto_mircuts_on 
    encoder_params.MIRCUTS = options.mircuts 
    encoder_params.AUTO_MODKCUTS_ON = options.auto_modkcuts_on 
    encoder_params.MODKCUTS = options.modkcuts 
    encoder_params.AUTO_NETWORKCUTS_ON = options.auto_networkcuts_on 
    encoder_params.NETWORKCUTS = options.networkcuts 
    encoder_params.AUTO_PROJIMPLIEDCUTS_ON = options.auto_projimpliedcuts_on 
    encoder_params.PROJIMPLIEDCUTS = options.projimpliedcuts 
    encoder_params.AUTO_RELAXLIFTCUTS_ON = options.auto_relaxliftcuts_on 
    encoder_params.RELAXLIFTCUTS = options.relaxliftcuts 
    encoder_params.AUTO_RLTCUTS_ON = options.auto_rltcuts_on 
    encoder_params.RLTCUTS = options.rltcuts 
    encoder_params.AUTO_STRONGCGCUTS_ON = options.auto_strongcgcuts_on 
    encoder_params.STRONGCGCUTS = options.strongcgcuts 
    encoder_params.AUTO_SUBMIPCUTS_ON = options.auto_submipcuts_on 
    encoder_params.SUBMIPCUTS = options.submipcuts 
    encoder_params.AUTO_ZEROHALFCUTS_ON = options.auto_zerohalfcuts_on 
    encoder_params.ZEROHALFCUTS = options.zerohalfcuts 
    encoder_params.AUTO_CUTPASSES_ON = options.auto_cutpasses_on
    encoder_params.CUTPASSES = options.cutpasses  

    encoder_params.XLAYER_DEP_CONSTRS = options.offline_dep
    encoder_params.LAYER_DEP_CONSTRS = options.offline_dep
    encoder_params.XLAYER_DEP_CUTS = options.online_dep
    encoder_params.LAYER_DEP_CUTS = options.online_dep
    encoder_params.IDEAL_CUTS = options.ideal_cuts
    encoder_params.WORKERS_NUMBER = options.workers

    # set splitter's parameters
    splitting_params = SplittingParam()
    splitting_params.FIXED_RATIO_CUTOFF = options.st_ratio
    splitting_params.DEPTH_POWER = options.depth_power
    splitting_params.STARTING_DEPTH_PARALLEL_SPLITTING_PROCESSES = options.splitters

    # create verifier
    verifier = VenusVerifier(options.net, spec, encoder_params, splitting_params, options.print)

    start = timer()
    result, job_id, extra = verifier.verify()
    end = timer()
    runtime = end - start

    result = verifier_result_to_string[result]

    print("Local robustness for input {} perturbed by {} over {}".format(options.lrob_input, options.lrob_radius, options.net),
          "is", result, "in {:9.4f}s".format(runtime), "job n", job_id)
    if result == True:
        ctx = extra.reshape(1, -1)
        network_output = nmodel.predict(x=ctx, batch_size=1)
        print("\t\tCounter-example:", list(extra))
        print("\t\tNetwork output :", list(network_output[0]))
        print("\t\tExpected label :", label)
    print("")

    return result, runtime, job_id


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

def main():
    parser = argparse.ArgumentParser(description="Venus Example",
                                     epilog="Exactly one of the parameters --acas_prop or --lrob_input is required.")
    parser.add_argument("--property", choices=["acas", "lrob"], default="lrob",
                        help="Verification property, one of acas or lrob (local robustness).")
    parser.add_argument("--net", type=str, default="/your/path/here/networks/RSL18a-linf01.h5",
                        help="Path to the neural network in Keras format.")
    parser.add_argument("--acas_prop", type=int, default=None,
                        help="Acas property number from 0 to 10. Default value is 1.")
    parser.add_argument("--lrob_input", type=str, default=None,
                        help="Path to the original input and the correct label for the local robustness property in the pickle format.")
    parser.add_argument("--lrob_radius", default=0.1, type=float,
                        help="Perturbation radius for L_inifinity norm. Default value is 0.1.")
    parser.add_argument("--st_ratio", default=0.5, type=float,
                        help="Cutoff value of the stable ratio during the splitting procedure. Default value is 0.5.")
    parser.add_argument("--depth_power", default=1.0, type=float,
                        help="Parameter for the splitting depth. Higher values favour splitting. Default value is 1.")
    parser.add_argument("--splitters", default=0, type=int,
                        help="Determines the number of splitting processes = 2^splitters. Default value is 0.")
    parser.add_argument("--workers", default=1, type=int,
                        help="Number of worker processes. Default value is 1.")
    parser.add_argument("--offline_dep", default=False, type=boolean_string,
                        help="Whether to include offline dependency cuts (before starting the solver) or not. Default value is True.")
    parser.add_argument("--online_dep", default=False, type=boolean_string,
                        help="Whether to include online dependency cuts (through solver callbacks) or not. Default value is True.")
    parser.add_argument("--ideal_cuts", default=False, type=boolean_string,
                        help="Whether to include online ideal cuts (through solver callbacks) or not. Default value is True.")
    parser.add_argument("--opt_bounds", default=False, type=boolean_string,
                        help="Whether to optimise bounds using linear relaxation or not. Default value is False.")
    parser.add_argument("--timeout", default=100000, type=int,
                        help="Timeout in seconds. Default value is 3600.")
    parser.add_argument("--logfile", default=None, type=str,
                        help="Path to logging file.")
    parser.add_argument("--print", default=False, type=boolean_string,
                        help="Print extra information or not. Default value is False.")

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

    ARGS = parser.parse_args()

    if ARGS.acas_prop is None and ARGS.lrob_input is None:
        parser.print_help()
        sys.exit(1)

    if ARGS.property == "acas":
        result, runtime, job_id = verify_acas_property(ARGS)
    else:
        result, runtime, job_id = verify_local_robustness(ARGS)

    if not ARGS.logfile is None:
        log = open(ARGS.logfile, 'a')
        log.write('{},'.format(ARGS.net.split(os.path.sep)[-1]))
        if ARGS.property == "acas":
            log.write('{},'.format(acas_properties[ARGS.acas_prop]['name']))
        else:
            log.write('{},{},'.format(ARGS.lrob_input.split(os.path.sep)[-1].split('.')[0], ARGS.lrob_radius))
        log.write('{},{:9.4f}\n'.format(result, runtime))
        log.close()

if __name__ == "__main__":
    main()
