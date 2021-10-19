from src.parameters import *
from src.layers import *
from gurobipy import *
from timeit import default_timer as timer
import numpy as np
from src.fldependencies import NegLiteral, FirstLayerDependencies, PosLiteral


def solver(vmodel):
    def _ineqs(model, layer, p_layer, node):
        # derive set of inequality nodes

        (s, e) = vmodel.get_var_indices(p_layer.depth, 'out')
        in_vars = model._vars[s:e]
        _in = np.asarray(model.cbGetNodeRel(in_vars))
        (s, e) = vmodel.get_var_indices(layer.depth, 'delta')
        delta_vars = model._vars[s:e]
        _delta = np.asarray(model.cbGetNodeRel(delta_vars))
        # in_vars = p_layer.vars['out'].tolist() 
        # _in = np.asarray(model.cbGetNodeRel(in_vars))
        # delta_vars = layer.vars['delta'].tolist()
        # _delta = np.asarray(model.cbGetNodeRel(delta_vars))
        ineqs = []

        for p_node in layer.getActiveWeights(node):
            l = layer.getPrimeLBound(node, p_layer, p_node)
            u = layer.getPrimeUBound(node, p_layer, p_node)
            lhs = layer.weights[node][p_node] * _in[p_node]
            rhs = layer.weights[node][p_node] * \
                  (l * (1 - _delta[node]) + u * _delta[node])
            if lhs < rhs:
                ineqs.append(p_node)
        return ineqs

    def _cond(model, ineqs, layer, p_layer, node):
        # check inequality condition  on ineqs 
        (s, e) = vmodel.get_var_indices(p_layer.depth, 'out')
        in_vars = model._vars[s:e]
        _in = np.asarray(model.cbGetNodeRel(in_vars))
        (s, e) = vmodel.get_var_indices(layer.depth, 'delta')
        delta_vars = model._vars[s:e]
        _delta = np.asarray(model.cbGetNodeRel(delta_vars))
        (s, e) = vmodel.get_var_indices(layer.depth, 'out')
        out_vars = model._vars[s:e]
        _out = np.asarray(model.cbGetNodeRel(out_vars))
        # in_vars = p_layer.vars['out'].tolist() 
        # _in = np.asarray(model.cbGetNodeRel(in_vars))
        # delta_vars = layer.vars['delta'].tolist()
        # _delta = np.asarray(model.cbGetNodeRel(delta_vars))
        # out_vars = layer.vars['out'].tolist()
        # _out = np.asarray(model.cbGetNodeRel(out_vars))
        s1 = 0
        s2 = 0

        for p_node in range(p_layer.output_shape):
            l = layer.getPrimeLBound(node, p_layer, p_node)
            u = layer.getPrimeUBound(node, p_layer, p_node)
            if p_node in ineqs:
                s1 += layer.weights[node][p_node] * \
                      (_in[p_node] - l * (1 - _delta[node]))
            else:
                s2 += layer.weights[node][p_node] * u * _delta[node]
        p = layer.bias[node] * _delta[node]
        if _out[node] > p + s1 + s2:
            return True
        else:
            return False

    def _constr(model, ineqs, layer, p_layer, node):
        # build constraint on ineqs

        (s, e) = vmodel.get_var_indices(p_layer.depth, 'out')
        in_vars = model._vars[s:e]
        (s, e) = vmodel.get_var_indices(layer.depth, 'delta')
        delta_vars = model._vars[s:e]
        (s, e) = vmodel.get_var_indices(layer.depth, 'out')
        out_vars = model._vars[s:e]
        # in_vars = p_layer.vars['out'].tolist() 
        # delta_vars = layer.vars['delta'].tolist()
        # out_vars = layer.vars['out'].tolist()
        le = LinExpr()
        s = 0

        for p_node in range(p_layer.output_shape):
            l = layer.getPrimeLBound(node, p_layer, p_node)
            u = layer.getPrimeUBound(node, p_layer, p_node)
            if p_node in ineqs:
                le.addTerms(layer.weights[node][p_node], in_vars[p_node])
                le.addConstant(- l * layer.weights[node][p_node])
                le.addTerms(l * layer.weights[node][p_node], delta_vars[node])
            else:
                s += layer.weights[node][p_node] * u
        le.addTerms(s + layer.bias[node], delta_vars[node])
        # le.addConstant(0.00001)

        return (out_vars[node], le)

    def _get_current_delta(model):
        delta = []
        _delta = []
        for i in vmodel.lmodel.layers:
            (s, e) = vmodel.get_var_indices(i.depth, 'delta')
            d = model._vars[s:e]
            _d = np.asarray(model.cbGetNodeRel(d))
            delta.append(d)
            _delta.append(_d)

        return [delta, _delta]

    def _get_lin_descr(delta, _delta):
        le = LinExpr()
        for i in range(len(vmodel.lmodel.layers)):
            d = delta[i]
            _d = _delta[i]
            for j in range(len(d)):
                if _d[j] == 0 and not vmodel.lmodel.layers[i].is_fixed(j):
                    le.addTerms(1, d[j])
                elif _d[j] == 1 and not vmodel.lmodel.layers[i].is_fixed(j):
                    le.addConstant(1)
                    le.addTerms(-1, d[j])
        return le

    def _freq(cut_type, depth=None):
        freq = 10000
        if vmodel.params.CALLBACK_FREQ == CallbackFreq.DEFAULT:
            if cut_type == Cuts.XLAYER:
                if vmodel.params.XLAYER_DEP_FREQ == CallbackFreq.LOG:
                    freq = math.ceil(math.log(GRB.Callback.MIPNODE_NODCNT + 1))
                elif vmodel.params.XLAYER_DEP_FREQ == CallbackFreq.POW:
                    freq = math.ceil(pow(GRB.Callback.MIPNODE_NODCNT + \
                                    1,vmodel.params.XLAYER_DEP_FREQ_CONST))
                elif vmodel.params.XLAYER_DEP_FREQ == CallbackFreq.CONST:
                    freq = vmodel.params.XLAYER_DEP_FREQ_CONST
            elif cut_type == Cuts.LAYER:
                if vmodel.params.LAYER_DEP_FREQ == CallbackFreq.LOG:
                    freq = math.ceil(math.log(GRB.Callback.MIPNODE_NODCNT + 1))
                elif vmodel.params.LAYER_DEP_FREQ == CallbackFreq.POW:
                    freq = math.ceil(pow(GRB.Callback.MIPNODE_NODCNT + \
                                         1,vmodel.params.LAYER_DEP_FREQ_CONST))
                elif vmodel.params.LAYER_DEP_FREQ == CallbackFreq.CONST:
                    freq = vmodel.params.LAYER_DEP_FREQ_CONST
            elif cut_type == Cuts.GROUP:
                if vmodel.params.GROUP_DEP_FREQ == CallbackFreq.LOG:
                    freq = math.ceil(math.log(GRB.Callback.MIPNODE_NODCNT + 1))
                elif vmodel.params.GROUP_DEP_FREQ == CallbackFreq.POW:
                    freq = math.ceil(pow(GRB.Callback.MIPNODE_NODCNT + \
                                         1, vmodel.params.GROUP_DEP_FREQ_CONST))
                elif vmodel.params.GROUP_DEP_FREQ == CallbackFreq.CONST:
                    freq = vmodel.params.GROUP_DEP_FREQ_CONST
            elif cut_type == Cuts.IDEAL:
                if vmodel.params.IDEAL_FREQ == CallbackFreq.LOG:
                    freq = math.ceil(math.log(GRB.Callback.MIPNODE_NODCNT + 1))
                elif vmodel.params.IDEAL_FREQ == CallbackFreq.POW:
                    freq = math.ceil(pow(GRB.Callback.MIPNODE_NODCNT + \
                                         1, vmodel.params.IDEAL_FREQ_CONST))
                elif vmodel.params.IDEAL_FREQ == CallbackFreq.CONST:
                    freq = vmodel.params.IDEAL_FREQ_CONST
        else:
            if vmodel.params.CALLBACK_FREQ == CallbackFreq.LOG:
                freq = math.ceil(math.log(GRB.Callback.MIPNODE_NODCNT + 1))
            elif vmodel.params.CALLBACK_FREQ == CallbackFreq.POW:
                freq = math.ceil(math.pow(GRB.Callback.MIPNODE_NODCNT \
                                          + 1, vmodel.params.CALLBACK_FREQ_CONST))
            elif vmodel.params.CALLBACK_FREQ == CallbackFreq.CONST:
                freq = vmodel.params.CALLBACK_FREQ_CONST

        if isinstance(depth, int):
            freq *= depth
        f = np.random.randint(0, freq, 1)
        if f == 0:
            return True
        else:
            return False

    def _compute_rnt_bounds(model):
        [delta, _delta] = _get_current_delta(model)
        le = _get_lin_descr(delta, _delta)
        if vmodel.params.DEP_DEPTH==-1:
            end = len(vmodel.lmodel.layers) - 1
        else:
            end = vmodel.params.DEP_DEPTH
        vmodel.compute_bounds(runtime=True, binary_vars=_delta, end=end)

        return le, delta, _delta

    def _add_dep_cuts(model):
        # calculate bounds and deps
        if vmodel.params.XLAYER_DEP_CUTS or \
        vmodel.params.LAYER_DEP_CUTS or vmodel.params.GROUP_DEP_CUTS:
            [delta, _delta] = _get_current_delta(model)
            le = _get_lin_descr(delta, _delta)
            if vmodel.params.DEP_DEPTH==-1:
                end = len(vmodel.lmodel.layers) - 1
            else:
                end = vmodel.params.DEP_DEPTH
            vmodel.compute_bounds(runtime=True, binary_vars=_delta, end=end)


            if vmodel.params.XLAYER_DEP_CUTS:
                _add_xlayer_dep_cuts(model, le, delta, _delta, end=end)
            if vmodel.params.LAYER_DEP_CUTS:
                _add_layer_dep_cuts(model, le, delta, _delta, end=end)
            if vmodel.params.GROUP_DEP_CUTS:
                _add_group_dep_cuts(model, le, delta, _delta, end=end)

    def _add_xlayer_dep_cuts(model, end=-1):

        # Add cuts as per the frequency parameter 
        if not _freq(Cuts.XLAYER,1):
            return
        # compute runtime bounds
        le, delta, _delta = _compute_rnt_bounds(model)
        # compute xlayer dependencies
        vmodel.compute_xlayer_deps(runtime=True, binary_vars=_delta, \
                                   end=end)
        # compute xlayer dependency cuts
        dep_cuts = 0
        for i in range(end - 1):
            l = vmodel.lmodel.layers[i]
            n_l = vmodel.lmodel.layers[i + 1]
            d = delta[i]
            n_d = delta[i + 1]
            if not isinstance(l,Relu) or not isinstance(n_l,Relu):
                continue
            ts = timer()
            for nd in range(l.output_shape):
                for (n_nd, dep) in l.xlayer_deps[nd]:
                    if dep == DepType.I_A:
                        model.cbCut(1 - n_d[n_nd] <= le + d[nd])
                    elif dep == DepType.I_I:
                        model.cbCut(n_d[n_nd] <= le + d[nd])
                    else:
                        raise Exception("""Unknown dependency type""", dep)
                    dep_cuts += 1

            te = timer()
            if vmodel.params.DEBUG_MODE:
                print(f'         - Added xlayer dependency cuts, #cuts: {dep_cuts}, layer: {l.depth}, callback time: {te - ts}')

    def _add_layer_dep_cuts(model, end=-1):

        # Add cuts as per the frequency parameter 
        if not _freq(Cuts.LAYER,1):
            return
        # compute runtime bounds
        le, delta, _delta = _compute_rnt_bounds(model)
        vmodel.compute_layer_deps(runtime=True, binary_vars=_delta, end=end)

        dep_cuts = 0
        p_l = vmodel.lmodel.input
        for i in range(end):
            l = vmodel.lmodel.layers[i]
            if not isinstance(l,Relu):
                p_l = l
                continue
            ts = timer()
            d = delta[i]
            for nd, (neg, pos) in l.layer_deps.dep_per_var.items():
                for x in neg:
                    if isinstance(x, PosLiteral):
                        model.cbCut(1 - d[x.i] <= le + d[nd])
                    else:
                        model.cbCut(d[x.i] <= le + d[nd])
                    dep_cuts += 1
                for x in pos:
                    if isinstance(x, PosLiteral):
                        model.cbCut(1 - d[x.i] <= le + 1 - d[nd])
                    else:
                        model.cbCut(d[x.i] <= le + 1 - d[nd])
                    dep_cuts += 1
            p_l = l
            te = timer()
            if vmodel.params.DEBUG_MODE:
                print(f'         - Added layer dependency cuts, #cuts: {dep_cuts}, layer: {l.depth}, callback time: {te - ts}')

    def _add_group_dep_cuts(model, end=-1):
        if not _freq(Cuts.GROUP, 1):
            return
        # compute group dependencies 
        vmodel.compute_group_deps(runtime=True, binary_vars=_delta, end=end)

        # add group dependency cuts
        dep_cuts = 0
        for i in range(end - 1):
            l = vmodel.lmodel.layers[i]
            n_l = vmodel.lmodel.layers[i + 1]
            d = delta[i]
            n_d = delta[i + 1]
            if not isinstance(l, Relu) or not isinstance(n_l, Relu):
                continue
            ts = timer()
            for nd in l.group_deps['group']:
                le.addTerms(1, d[nd])
            for (n_nd, dep) in l.group_deps['nodes']:
                if dep == DepType.I_A:
                    model.cbCut(1 - n_d[n_nd] <= le)
                elif dep == DepType.I_I:
                    model.cbCut(n_d[n_nd] <= le)
                else:
                    raise Exception("""Unknown dependency type""")
                dep_cuts += 1
            te = timer()
            if vmodel.params.DEBUG_MODE:
                print(f'        - Added group dependency cuts, #cuts: {dep_cuts}, layer: {l.depth}, callback time: {te - ts}')

    def _add_dep_cuts_old(model):
        ts = timer()
        # calculate bounds and deps
        delta = []
        le = LinExpr()
        for i in vmodel.lmodel.layers:
            (s, e) = vmodel.get_var_indices(i.depth, 'delta')
            d = model.vars[s:e]
            _d = np.asarray(model.cbGetNodeRel(d))
            for j in range(len(d)):
                if _d[j] == 0 and not i.is_fixed(j):
                    le.addTerms(1, d[j])
                elif _d[j] == 1 and not i.is_fixed(j):
                    le.addConstant(1)
                    le.addTerms(-1, d[j])
            delta.append(_d)
        vmodel.compute_bounds(runtime=True, binary_vars=delta, end=vmodel.params.DEP_DEPTH)
        vmodel.compute_xlayer_deps(runtime=True, binrary_vars=delta)

        # add dependency cuts
        dep_cuts = 0
        for i in range(len(vmodel.lmodel.layers) - 1):
            l = vmodel.lmodel.layers[i]
            n_l = vmodel.lmodel.layers[i + 1]
            (s, e) = vmodel.get_var_indices(l.depth, 'delta')
            delta = model.vars[s:e]
            (s, e) = vmodel.get_var_indices(n_l.depth, 'delta')
            n_delta = model.vars[s:e]
            if isinstance(l,Relu):
                for nd in range(l.output_shape):
                    for (n_nd, dep) in l.deps[nd]:
                        if dep == DepType.I_A:
                            model.cbCut(1 - n_delta[n_nd] <= le + delta[nd])
                        elif dep == DepType.I_I:
                            model.cbCut(n_delta[n_nd] <= le + delta[nd])
                        else:
                            print('\n\n here \n\n')
                        dep_cuts += 1
        te = timer()
        if vmodel.params.DEBUG_MODE:
            print(f'         - Added dependency cuts, #cuts: {dep_cuts}, callback time: {te - ts}')

    def _add_ideal_cuts(model):
        p_l = vmodel.lmodel.input
        for l in vmodel.lmodel.layers:
            if not isinstance(l,ReluIdeal):
                p_l = l
                continue
            if not _freq(Cuts.IDEAL,l.depth):
                p_l = l
                continue 
            ts = timer()
            for nd in range(l.output_shape):
                # only check not fixed nodes
                if not l.is_fixed(nd):
                    # only check nodes where 0 < delta < 1 
                    (s, e) = vmodel.get_var_indices(l.depth, 'delta')
                    delta = model._vars[s:e]
                    _delta = model.cbGetNodeRel(delta)
                    if _delta[nd] > 0 and _delta[nd] < 1:
                        ineqs = _ineqs(model, l, p_l, nd)
                        if _cond(model, ineqs, l, p_l, nd):
                            (lhs, rhs) = _constr(model, ineqs, l, p_l, nd)
                            model.cbCut(lhs <= rhs)
            te = timer()
            if vmodel.params.DEBUG_MODE:
                print(f'         - Added ideal cuts, layer: {l.depth}, callback time: { te - ts}')
            p_l = l

    def _callback(model, where):
        if where == GRB.Callback.MIPNODE:
            # terminate if a robustness cex has already been found
            if vmodel.spec.isMaxLRob() and \
                    model.cbGet(GRB.Callback.MIPNODE_OBJBST) > 0:
                model.terminate()
            if model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.Status.OPTIMAL:
                if vmodel.params.DEP_DEPTH==-1:
                    end = len(vmodel.lmodel.layers) - 1
                else:
                    end = vmodel.params.DEP_DEPTH
                if vmodel.params.IDEAL_CUTS:
                    _add_ideal_cuts(model)
                if vmodel.params.XLAYER_DEP_CUTS:
                    _add_xlayer_dep_cuts(model,end)
                if vmodel.params.LAYER_DEP_CUTS:
                    _add_layer_dep_cuts(model,end)
                if vmodel.params.GROUP_DEP_CUTS:
                    _add_group_dep_cuts(model,end)

              
    # Set gurobi parameters
    if not vmodel.params.TIME_LIMIT == -1:
        vmodel.gmodel.setParam('TIME_LIMIT', vmodel.params.TIME_LIMIT)
    vmodel.gmodel.setParam('aggfill', vmodel.params.AGGFILL)
    vmodel.gmodel.setParam('aggregate', vmodel.params.AGGREGATE)
    vmodel.gmodel.setParam('auto_aggfill_on', vmodel.params.AUTO_AGGFILL_ON)
    vmodel.gmodel.setParam('auto_barcorrectors', vmodel.params.AUTO_BARCORRECTORS)
    vmodel.gmodel.setParam('auto_barhomogeneous', vmodel.params.AUTO_BARHOMOGENEOUS)
    vmodel.gmodel.setParam('auto_barorder', vmodel.params.AUTO_BARORDER)
    vmodel.gmodel.setParam('auto_bqpcuts_on', vmodel.params.AUTO_BQPCUTS_ON)
    vmodel.gmodel.setParam('auto_branchdir_on', vmodel.params.AUTO_BRANCHDIR_ON)
    vmodel.gmodel.setParam('auto_cliquecuts_on', vmodel.params.AUTO_CLIQUECUTS_ON)
    vmodel.gmodel.setParam('auto_covercuts_on', vmodel.params.AUTO_COVERCUTS_ON)
    vmodel.gmodel.setParam('auto_crossover', vmodel.params.AUTO_CROSSOVER)
    vmodel.gmodel.setParam('auto_cutaggpasses_on', vmodel.params.AUTO_CUTAGGPASSES_ON)
    vmodel.gmodel.setParam('auto_cutpasses_on', vmodel.params.AUTO_CUTPASSES_ON)
    vmodel.gmodel.setParam('auto_cuts_on', vmodel.params.AUTO_CUTS_ON)
    vmodel.gmodel.setParam('auto_degenmoves_on', vmodel.params.AUTO_DEGENMOVES_ON)
    vmodel.gmodel.setParam('auto_disconnected_on', vmodel.params.AUTO_DISCONNECTED_ON)
    vmodel.gmodel.setParam('auto_flowcovercuts_on', vmodel.params.AUTO_FLOWCOVERCUTS_ON)
    vmodel.gmodel.setParam('auto_flowpathcuts_on', vmodel.params.AUTO_FLOWPATHCUTS_ON)
    vmodel.gmodel.setParam('auto_gomorypasses_on', vmodel.params.AUTO_GOMORYPASSES_ON)
    vmodel.gmodel.setParam('auto_gubcovercuts_on', vmodel.params.AUTO_GUBCOVERCUTS_ON)
    vmodel.gmodel.setParam('auto_impliedcuts_on', vmodel.params.AUTO_IMPLIEDCUTS_ON)
    vmodel.gmodel.setParam('auto_infproofcuts_on', vmodel.params.AUTO_INFPROOFCUTS_ON)
    vmodel.gmodel.setParam('auto_issmethod_on', vmodel.params.AUTO_ISSMETHOD_ON)
    vmodel.gmodel.setParam('auto_minrelnodes_on', vmodel.params.AUTO_MINRELNODES_ON)
    vmodel.gmodel.setParam('auto_mipsepcuts_on', vmodel.params.AUTO_MIPSEPCUTS_ON)
    vmodel.gmodel.setParam('auto_mircuts_on', vmodel.params.AUTO_MIRCUTS_ON)
    vmodel.gmodel.setParam('auto_modkcuts_on', vmodel.params.AUTO_MODKCUTS_ON)
    vmodel.gmodel.setParam('auto_networkcuts_on', vmodel.params.AUTO_NETWORKCUTS_ON)
    vmodel.gmodel.setParam('auto_nodemethod_on', vmodel.params.AUTO_NODEMETHOD_ON)
    vmodel.gmodel.setParam('auto_normadjust_on', vmodel.params.AUTO_NORMADJUST_ON)
    vmodel.gmodel.setParam('auto_predual_on', vmodel.params.AUTO_PREDUAL_ON)
    vmodel.gmodel.setParam('auto_prepasses_on', vmodel.params.AUTO_PREPASSES_ON)
    vmodel.gmodel.setParam('auto_presolve_on', vmodel.params.AUTO_PRESOLVE_ON)
    vmodel.gmodel.setParam('auto_presparsify_on', vmodel.params.AUTO_PRESPARSIFY_ON)
    vmodel.gmodel.setParam('auto_projimpliedcuts_on', vmodel.params.AUTO_PROJIMPLIEDCUTS_ON)
    vmodel.gmodel.setParam('auto_pumppasses_on', vmodel.params.AUTO_PUMPPASSES_ON)
    vmodel.gmodel.setParam('auto_quad_on', vmodel.params.AUTO_QUAD_ON)
    vmodel.gmodel.setParam('auto_relaxliftcuts_on', vmodel.params.AUTO_RELAXLIFTCUTS_ON)
    vmodel.gmodel.setParam('auto_rins_on', vmodel.params.AUTO_RINS_ON)
    vmodel.gmodel.setParam('auto_rltcuts_on', vmodel.params.AUTO_RLTCUTS_ON)
    vmodel.gmodel.setParam('auto_sifting_on', vmodel.params.AUTO_SIFTING_ON)
    vmodel.gmodel.setParam('auto_siftmethod_on', vmodel.params.AUTO_SIFTMETHOD_ON)
    vmodel.gmodel.setParam('auto_simplexpricing_on', vmodel.params.AUTO_SIMPLEXPRICING_ON)
    vmodel.gmodel.setParam('auto_startnodelimit', vmodel.params.AUTO_STARTNODELIMIT)
    vmodel.gmodel.setParam('auto_strongcgcuts_on', vmodel.params.AUTO_STRONGCGCUTS_ON)
    vmodel.gmodel.setParam('auto_submipcuts_on', vmodel.params.AUTO_SUBMIPCUTS_ON)
    vmodel.gmodel.setParam('auto_symmetry_on', vmodel.params.AUTO_SYMMETRY_ON)
    vmodel.gmodel.setParam('auto_varbranch_on', vmodel.params.AUTO_VARBRANCH_ON)
    vmodel.gmodel.setParam('auto_zerohalfcuts_on', vmodel.params.AUTO_ZEROHALFCUTS_ON)
    vmodel.gmodel.setParam('auto_zeroobjnodes_on', vmodel.params.AUTO_ZEROOBJNODES_ON)
    vmodel.gmodel.setParam('barcorrectors', vmodel.params.BARCORRECTORS)
    vmodel.gmodel.setParam('barorder', vmodel.params.BARORDER)
    vmodel.gmodel.setParam('crossoverbasis', vmodel.params.CROSSOVERBASIS)
    vmodel.gmodel.setParam('cuts', vmodel.params.CUTS)
    vmodel.gmodel.setParam('disconnected', vmodel.params.DISCONNECTED)
    vmodel.gmodel.setParam('dualreductions', vmodel.params.DUALREDUCTIONS)
    vmodel.gmodel.setParam('heuristics', vmodel.params.HEURISTICS)
    vmodel.gmodel.setParam('improvestartgap', vmodel.params.IMPROVESTARTGAP)
    vmodel.gmodel.setParam('improvestartnodes', vmodel.params.IMPROVESTARTNODES)
    vmodel.gmodel.setParam('improvestarttime', vmodel.params.IMPROVESTARTTIME)
    vmodel.gmodel.setParam('infproofcuts', vmodel.params.INFPROOFCUTS)
    vmodel.gmodel.setParam('infunbdinfo', vmodel.params.INFUNBDINFO)
    vmodel.gmodel.setParam('mipfocus', vmodel.params.MIPFOCUS)
    vmodel.gmodel.setParam('networkcuts', vmodel.params.NETWORKCUTS)
    vmodel.gmodel.setParam('partitionplace', vmodel.params.PARTITIONPLACE)
    vmodel.gmodel.setParam('perturbvalue', vmodel.params.PERTURBVALUE)
    vmodel.gmodel.setParam('precrush', vmodel.params.PRECRUSH)
    vmodel.gmodel.setParam('predeprow', vmodel.params.PREDEPROW)
    vmodel.gmodel.setParam('prepasses', vmodel.params.PREPASSES)
    vmodel.gmodel.setParam('presparsify', vmodel.params.PRESPARSIFY)
    vmodel.gmodel.setParam('pumppasses', vmodel.params.PUMPPASSES)
    vmodel.gmodel.setParam('quad', vmodel.params.QUAD)
    vmodel.gmodel.setParam('siftmethod', vmodel.params.SIFTMETHOD)
    vmodel.gmodel.setParam('simplexpricing', vmodel.params.SIMPLEXPRICING)
    vmodel.gmodel.setParam('startnodelimit', vmodel.params.STARTNODELIMIT)
    vmodel.gmodel.setParam('submipcuts', vmodel.params.SUBMIPCUTS)
    vmodel.gmodel.setParam('submipnodes', vmodel.params.SUBMIPNODES)
    vmodel.gmodel.setParam('symmetry', vmodel.params.SYMMETRY)
    vmodel.gmodel.setParam('threads', vmodel.params.THREADS)

    vmodel.gmodel._vars = vmodel.gmodel.getVars()
    if not vmodel.params.DEFAULT_CUTS:
        vmodel.disable_cuts()

    # Optimise
    if vmodel.params.IDEAL_CUTS or vmodel.params.XLAYER_DEP_CUTS or \
    vmodel.params.LAYER_DEP_CUTS or vmodel.params.GROUP_DEP_CUTS:
        vmodel.gmodel.optimize(_callback)
    else:
        vmodel.gmodel.optimize()

    if vmodel.gmodel.status == GRB.OPTIMAL:
        cex_shape = vmodel.lmodel.input.vars['out'].shape
        cex = np.empty(shape=cex_shape)
        for i in itertools.product(*[range(j) for j in cex_shape]):
            cex[i] = vmodel.lmodel.input.vars['out'][i].x

        # import matplotlib.pyplot as plt
        # imgplot = plt.imshow(cex.reshape(28,28))
        # plt.show()

        return (True, cex)
        # return (True, vmodel.gmodel.MIPGap)
    elif vmodel.gmodel.status == GRB.TIME_LIMIT:
        return ('Timeout', None)
    elif vmodel.gmodel.status == GRB.INTERRUPTED:
        return ('Interrupted', None)
    elif vmodel.gmodel.status == GRB.INFEASIBLE or \
            vmodel.gmodel.status == GRB.INF_OR_UNBD:
        # vmodel.gmodel.computeIIS()
        # vmodel.gmodel.write('incons.ilp')
        return (False, None)

    return ("Unknown status {}".format(vmodel.gmodel.status), None)
