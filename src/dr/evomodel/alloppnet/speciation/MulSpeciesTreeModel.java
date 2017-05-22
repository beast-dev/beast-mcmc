/*
 * MulSpeciesTreeModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.alloppnet.speciation;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.evomodel.alloppnet.operators.MulTreeNodeSlide;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.alloppnet.parsers.MulSpeciesTreeModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.Scalable;
import dr.evomodel.alloppnet.util.AlloppMisc;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.HeapSort;
import jebl.util.FixedBitSet;

import java.util.*;
import java.util.logging.Logger;

/**
 * Multiply labelled species tree which includes demographic function per branch.
 *
 * @author Joseph Heled, Graham Jones
 *         Date: 20/12/2011
 */

/*
 * grjtodo-oneday This includes a lot of code from SpeciesTreeModel. It would be better to re-use!
 * Not clear how: neither contain nor extend seems ideal. This uses
 * a MulSpeciesBindings which is seriously different from a SpeciesBindings.
 * 
 * Or, at least simplify this, removing the bmp = brownian motion prior
 * and the piece-wise linear model for population sizes. 
 * 
 */


public class MulSpeciesTreeModel extends AbstractModel implements MutableTree, Citable, TreeTraitProvider, TreeLogger.LogUpon, Scalable {
    private final SimpleTree spTree;
    private final MulSpeciesBindings mulspb;
    private final Map<NodeRef, NodeProperties> props = new HashMap<NodeRef, NodeProperties>();
    public final Parameter sppSplitPopulations;
    private int[] singleStartPoints;
    private int[] pairStartPoints;

    private final Parameter coalPointsPops;
    private final Parameter coalPointsIndicator;
    private boolean nodePropsReady;

    private final NodeRef[] children;
    private final double[] heights;

    // any change of underlying parameters / models
    private boolean anyChange;
    // Tree has been edited in this cycle
    private boolean treeChanged;

    private final String spIndexAttrName = "spi";
    private final boolean bmp;
    private final boolean nonConstRootPopulation;
    private final boolean constantPopulation;

    // grj
    public final static boolean DBUGTUNE = false;



    private class NodeProperties {
        private final int speciesIndex;
        public VDdemographicFunction demogf;
        FixedBitSet spSet;

        public NodeProperties(int n) {
            speciesIndex = n;
            demogf = null;
            spSet = new FixedBitSet(mulspb.nSpSeqs());
        }
    }

    public MulSpeciesTreeModel(MulSpeciesBindings species, Parameter sppSplitPopulations,
                               Parameter coalPointsPops, Parameter coalPointsIndicator, Tree startTree,
                               boolean bmp, boolean nonConstRootPopulation, boolean constantPopulation) {
        super(MulSpeciesTreeModelParser.MUL_SPECIES_TREE);

        this.mulspb = species;

        this.sppSplitPopulations = sppSplitPopulations;
        this.coalPointsPops = coalPointsPops;
        this.coalPointsIndicator = coalPointsIndicator;

        this.bmp = bmp;
        this.nonConstRootPopulation = nonConstRootPopulation;
        this.constantPopulation = constantPopulation;

        addVariable(sppSplitPopulations);

        addModel(species);

        if (coalPointsPops != null) {
            assert coalPointsIndicator != null;
            assert !constantPopulation;

            addVariable(coalPointsPops);
            addVariable(coalPointsIndicator);

            final double[][] pts = species.getPopTimesSingle();

            int start = 0;
            singleStartPoints = new int[pts.length];
            for (int i = 0; i < pts.length; i++) {
                singleStartPoints[i] = start;
                start += pts[i].length;
            }

            if (!bmp) {
                final double[][] ptp = species.getPopTimesPair();
                pairStartPoints = new int[ptp.length];
                for (int i = 0; i < ptp.length; i++) {
                    pairStartPoints[i] = start;
                    start += ptp[i].length;
                }
            }
        }

        // build an initial noninformative tree
        spTree = compatibleUninformedSpeciesTree(startTree);

        // some of the code is generic but some parts assume a binary tree.
        assert TreeUtils.isBinary(spTree);

        final int nNodes = spTree.getNodeCount();
        heights = new double[nNodes];
        children = new NodeRef[2 * nNodes + 1];

        // fixed properties
        for (int k = 0; k < getExternalNodeCount(); ++k) {
            final NodeRef nodeRef = getExternalNode(k);
            final int n = (Integer) getNodeAttribute(nodeRef, spIndexAttrName);
            final NodeProperties np = new NodeProperties(n);
            props.put(nodeRef, np);
            np.spSet.set(n);
        }

        for (int k = 0; k < getInternalNodeCount(); ++k) {
            final NodeRef nodeRef = getInternalNode(k);
            props.put(nodeRef, new NodeProperties(-1));
        }

        nodePropsReady = false;

        // crappy way to pass a result back from compatibleUninformedSpeciesTree.
        // check is using isCompatible(), which requires completion of construction.
        boolean check = spTree.getAttribute("check") != null;
        spTree.setAttribute("check", null);

        if (check) {
            // change only is needed - if the user provided a compatible state she may know
            // what she is doing
            for (MulSpeciesBindings.GeneTreeInfo gt : species.getGeneTrees()) {
                if (!isCompatible(gt)) {

                    species.makeCompatible(spTree.getRootHeight());
                    for (MulSpeciesBindings.GeneTreeInfo t : species.getGeneTrees()) {
                        assert isCompatible(t);
                    }
                    anyChange = false;
                    break;
                }
            }
        }

        Logger.getLogger("dr.evomodel.speciation.allopolyploid").info("\tConstructing a multiply labelled tree, please cite:\n"
                + Citable.Utils.getCitationString(this));

    }



    // grj
    private String nodeAsText(NodeRef node, int indentlen) {
        StringBuilder s = new StringBuilder();
        Formatter formatter = new Formatter(s, Locale.US);
        if (spTree.isExternal(node)) {
            formatter.format("%s ", spTree.getNodeTaxon(node));
        } else {
            formatter.format("%s ", "+");
        }
        while (s.length() < 20-indentlen) {
            formatter.format("%s", " ");
        }
        formatter.format("%s ", AlloppMisc.nonnegIn8Chars(spTree.getNodeHeight(node)));
        // it would be nice to disply popsizes and nlineages like allopp mul tree
        return s.toString();
    }



    // grj
    private String subtreeAsText(NodeRef node, String s, Stack<Integer> x, int depth, String b) {
        Integer[] y = x.toArray(new Integer[x.size()]);
        StringBuffer indent = new StringBuffer();
        for (int i = 0; i < depth; i++) {
            indent.append("  ");
        }
        for (int i = 0; i < y.length; i++) {
            indent.replace(2*y[i], 2*y[i]+1, "|");
        }
        if (b.length() > 0) {
            indent.replace(indent.length()-b.length(), indent.length(), b);
        }
        s += indent;
        s += nodeAsText(node, indent.length());
        s += System.getProperty("line.separator");
        String subs = "";
        if (!spTree.isExternal(node)) {
            x.push(depth);
            subs += subtreeAsText(spTree.getChild(node, 0), "", x, depth+1, "-");
            x.pop();
            subs += subtreeAsText(spTree.getChild(node, 1), "", x, depth+1, "`-");
        }
        return s + subs;
    }



    // grj
    public String asText() {
        String header = "topology             height" + System.getProperty("line.separator");

        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + subtreeAsText(spTree.getRoot(), s, x, 0, "");
    }


    public String toString() {
        int ngt = mulspb.numberOfGeneTrees();
        String nl = System.getProperty("line.separator");
        String s = nl + asText() + nl;
        for (int g = 0; g < ngt; g++) {
            s += "Gene tree " + g + nl;
            s += mulspb.genetreeAsText(g) + nl;
            s += mulspb.seqassignsAsText(g) + nl;
        }
        s += nl;
        return s;
    }



    // grj
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[1];
        columns[0] = new LogColumn.Default("    MUL-tree and gene trees", this);
        return columns;
    }





    public boolean constPopulation() {
        return constantPopulation;
    }

    // Is gene tree compatible with species tree
    // grj: the sequence assignments are handled in geneTreeInfo.getCoalInfo()
    public boolean isCompatible(MulSpeciesBindings.GeneTreeInfo geneTreeInfo) {
        // can't set demographics if a tree is not compatible, but we need spSets.
        if (!nodePropsReady) {
            setSPsets(getRoot());
        }
        return isSubtreeCompatible(getRoot(), geneTreeInfo.getCoalInfo(), 0) >= 0;
    }



    // Not very efficient, should do something better, based on traversing the cList once
    private int isSubtreeCompatible(NodeRef node, MulSpeciesBindings.CoalInfo[] cList, int loc) {
        if (!isExternal(node)) {
            int l = -1;
            for (int nc = 0; nc < getChildCount(node); ++nc) {
                int l1 = isSubtreeCompatible(getChild(node, nc), cList, loc);
                if (l1 < 0) {
                    return -1;
                }
                assert l == -1 || l1 == l;

                l = l1;
            }
            loc = l;

            assert cList[loc].ctime >= getNodeHeight(node);
        }

        if (node == getRoot()) {
            return cList.length;
        }

        // spSet guaranteed to be ready by caller
        final FixedBitSet nodeSps = props.get(node).spSet;

        final double limit = getNodeHeight(getParent(node));

        while (loc < cList.length) {
            final MulSpeciesBindings.CoalInfo ci = cList[loc];
            if (ci.ctime >= limit) {
                break;
            }
            boolean allIn = true, noneIn = true;

            for (int i = 0; i < 2; ++i) {
                final FixedBitSet s = ci.sinfo[i];

                final int in1 = s.intersectCardinality(nodeSps);
                if (in1 > 0) {
                    noneIn = false;
                }
                if (s.cardinality() != in1) {
                    allIn = false;
                }
            }
            if (!(allIn || noneIn)) {
                return -1;
            }
            ++loc;
        }
        return loc;
    }


    private static double
    fp(double val, double low, double[][] tt, int[] ii) {
        for (int k = 0; k < ii.length; ++k) {
            int ip = ii[k];
            if (ip == tt[k].length || val <= tt[k][ip]) {
                --ip;
                while (ip >= 0 && val <= tt[k][ip]) {
                    --ip;
                }
                assert ((ip < 0) || (tt[k][ip] < val)) && ((ip + 1 == tt[k].length) || (val <= tt[k][ip + 1]));
                if (ip >= 0) {
                    low = Math.max(low, tt[k][ip]);
                }
            } else {
                ++ip;
                while (ip < tt[k].length && val > tt[k][ip]) {
                    ++ip;
                }
                assert tt[k][ip - 1] < val && ((ip == tt[k].length) || (val <= tt[k][ip]));
                low = Math.max(low, tt[k][ip - 1]);
            }
        }
        return low;
    }

    private interface SimpleDemographicFunction {
        double population(double t);

        double upperBound();
    }

    private class PLSD implements SimpleDemographicFunction {
        private final double[] pops;
        private final double[] times;

        public PLSD(double[] pops, double[] times) {
            assert pops.length == times.length + 1;

            this.pops = pops;
            this.times = times;
        }

        public double population(double t) {
            if (t >= upperBound()) {
                return pops[pops.length - 1];
            }

            int k = 0;
            while (t > times[k]) {
                t -= times[k];
                ++k;
            }
            double a = t / (times[k] - (k > 0 ? times[k - 1] : 0));
            return a * pops[k] + (1 - a) * pops[k + 1];
        }

        public double upperBound() {
            return times[times.length - 1];
        }
    }
    // Pass arguments of recursive functions in a compact format.

    private class Args {
        final double[][] cps = mulspb.getPopTimesSingle();
        final double[][] cpp; // = species.getPopTimesPair();
        final int[] iSingle = new int[cps.length];
        final int[] iPair; // = new int[cpp.length];

        final double[] indicators = ((Parameter.Default) coalPointsIndicator).inspectParameterValues();
        final double[] pops = ((Parameter.Default) coalPointsPops).inspectParameterValues();
        final SimpleDemographicFunction[] dms;

        Args(Boolean bmp) {
            if (!bmp) {
                cpp = mulspb.getPopTimesPair();
                iPair = new int[cpp.length];
                dms = null;
            } else {
                cpp = null;
                iPair = null;
                int nsps = cps.length;
                dms = new SimpleDemographicFunction[nsps];
                for (int nsp = 0; nsp < nsps; ++nsp) {
                    final int start = singleStartPoints[nsp];
                    final int stop = nsp < nsps - 1 ? singleStartPoints[nsp + 1] : pops.length;

                    double[] pop = new double[1 + stop - start];
                    pop[0] = sppSplitPopulations.getParameterValue(nsp); //  pops[nsp];

                    for (int k = 0; k < stop - start; ++k) {
                        pop[k + 1] = pops[start + k];
                    }
                    dms[nsp] = new PLSD(pop, cps[nsp]);
                }
            }
        }

        private double findPrev(double val, double low) {
            low = fp(val, low, cps, iSingle);
            low = fp(val, low, cpp, iPair);

            return low;
        }
    }

    class RawPopulationHelper {
        final int[] preOrderIndices = new int[getNodeCount()];
        final double[] pops = ((Parameter.Default) sppSplitPopulations).inspectParameterValues();
        final int nsp = mulspb.nSpSeqs();
        final Args args = coalPointsPops != null ? new Args(bmp) : null;

        RawPopulationHelper() {
            setPreorderIndices(preOrderIndices);
        }

        public void getPopulations(NodeRef n, int nc, double[] p) {
            p[1] = pops[nsp + 2 * preOrderIndices[n.getNumber()] + nc];

            final NodeRef child = getChild(n, nc);
            if (isExternal(child)) {
                p[0] = pops[props.get(child).speciesIndex];
            } else {
                int k = nsp + 2 * preOrderIndices[child.getNumber()];
                p[0] = pops[k] + pops[k + 1];
            }
        }

        public double tipPopulation(NodeRef tip) {
            return pops[props.get(tip).speciesIndex];
        }

        // grj: this doesn't work for multiply labelled tree. 
        // Use nSqSeqs() instead for most purposes 
        /*public int nSpecies() {
            return species.nSpecies();
        }*/

        public boolean perSpeciesPopulation() {
            return args != null;
        }

        public double[] getTimes(int ns) {
            return ((PLSD) args.dms[ns]).times;
        }

        public double[] getPops(int ns) {
            return ((PLSD) args.dms[ns]).pops;
        }

        public void getRootPopulations(double[] p) {
            int k = nsp + 2 * preOrderIndices[getRoot().getNumber()];
            p[0] = pops[k] + pops[k + 1];
            p[1] = nonConstRootPopulation ? pops[pops.length - 1] : p[0];
        }

        public double geneTreesRootHeight() {
            //getNodeDemographic(getRoot()).
            double h = -1;
            for (MulSpeciesBindings.GeneTreeInfo t : mulspb.getGeneTrees()) {
                h = Math.max(h, t.tree.getNodeHeight(t.tree.getRoot()));
            }
            return h;
        }
    }

    RawPopulationHelper getPopulationHelper() {
        return new RawPopulationHelper();
    }

    static private class Points implements Comparable<Points> {
        final double time;
        double population;
        final boolean use;

        Points(double t, double p) {
            time = t;
            population = p;
            use = true;
        }

        Points(double t, boolean u) {
            time = t;
            population = 0;
            use = u;
        }

        public int compareTo(Points points) {
            return time < points.time ? -1 : (time > points.time ? 1 : 0);
        }
    }

    private NodeProperties
    setSPsets(NodeRef nodeID) {
        final NodeProperties nprop = props.get(nodeID);

        if (!isExternal(nodeID)) {
            nprop.spSet = new FixedBitSet(mulspb.nSpSeqs());
            for (int nc = 0; nc < getChildCount(nodeID); ++nc) {
                NodeProperties p = setSPsets(getChild(nodeID, nc));
                nprop.spSet.union(p.spSet);
            }
        }
        return nprop;
    }

    private int ti2f(int i, int j) {
        return (i == 0) ? j : 2 * i + j + 1;
    }

    private VDdemographicFunction
    bestLinearFit(double[] xs, double[] ys, boolean[] use) {

        assert (xs.length + 1) == ys.length;
        assert ys.length == use.length + 2 || ys.length == use.length + 1;

        int N = ys.length;
        if (N == 2) {
            // cheaper
            return new VDdemographicFunction(xs, ys, getUnits());
        }

        List<Integer> iv = new ArrayList<Integer>(2);
        iv.add(0);
        for (int k = 0; k < N - 2; ++k) {
            if (use[k]) {
                iv.add(k + 1);
            }
        }
        iv.add(N - 1);

        double[] ati = new double[xs.length + 1];
        ati[0] = 0.0;
        System.arraycopy(xs, 0, ati, 1, xs.length);
        int n = iv.size();

        double[] a = new double[3 * n];
        double[] v = new double[n];

        for (int k = 0; k < n - 1; ++k) {
            int i0 = iv.get(k);
            int i1 = iv.get(k + 1);

            double u0 = ati[i0];
            double u1 = ati[i1] - ati[i0];
            // on last interval add data for last point
            if (i1 == N - 1) {
                i1 += 1;
            }

            final int l = ti2f(k, k);
            final int l1 = ti2f(k + 1, k);

            for (int j = i0; j < i1; ++j) {
                double t = ati[j];
                double y = ys[j];

                double z = (t - u0) / u1;
                v[k] += y * (1 - z);

                a[l] += (1 - z) * (1 - z);
                a[l + 1] += z * (1 - z);

                a[l1] += z * (1 - z);
                a[l1 + 1] += z * z;
                v[k + 1] += y * z;
            }
        }

        for (int k = 0; k < n - 1; ++k) {

            final double r = a[ti2f(k + 1, k)] / a[ti2f(k, k)];
            for (int j = k; j < k + 3; ++j) {
                a[ti2f((k + 1), j)] -= a[ti2f(k, j)] * r;
            }
            v[k + 1] -= v[k] * r;
        }

        double[] z = new double[n];
        for (int k = n - 1; k > 0; --k) {
            z[k] = v[k] / a[ti2f(k, k)];
            v[k - 1] -= a[ti2f((k - 1), k)] * z[k];
        }

        z[0] = v[0] / a[ti2f(0, 0)];
        double[] t = new double[iv.size() - 1];
        for (int j = 0; j < t.length; ++j) {
            t[j] = ati[iv.get(j + 1)];
        }
        return new VDdemographicFunction(t, z, getUnits());
    }

    //  Assign positions in 'pointsList' for the sub-tree rooted at the ancestor of
    //  nodeID.
    //
    //  pointsList is indexed by node-id. Every element is a list of internal
    //  population points for the branch between nodeID and it's ancestor
    //

    private NodeProperties getDemographicPoints(final NodeRef nodeID, Args args, Points[][] pointsList) {

        final NodeProperties nprop = props.get(nodeID);
        final int nSpSeqs = mulspb.nSpSeqs();

        // Species assignment from the tips never changes
        if (!isExternal(nodeID)) {
            nprop.spSet = new FixedBitSet(nSpSeqs);
            for (int nc = 0; nc < getChildCount(nodeID); ++nc) {
                final NodeProperties p = getDemographicPoints(getChild(nodeID, nc), args, pointsList);
                nprop.spSet.union(p.spSet);
            }
        }

        if (args == null) {
            return nprop;
        }

        // parent height
        final double cHeight = nodeID != getRoot() ? getNodeHeight(getParent(nodeID)) : Double.MAX_VALUE;

        // points along branch
        // not sure what a good default size is?
        List<Points> allPoints = new ArrayList<Points>(5);

        if (bmp) {
            for (int isp = nprop.spSet.nextOnBit(0); isp >= 0; isp = nprop.spSet.nextOnBit(isp + 1)) {
                final double[] cp = args.cps[isp];
                final int upi = singleStartPoints[isp];

                int i = args.iSingle[isp];
                for (/**/; i < cp.length && cp[i] < cHeight; ++i) {
                    allPoints.add(new Points(cp[i], args.indicators[upi + i] > 0));
                }
                args.iSingle[isp] = i;
            }
        } else {

            for (int isp = nprop.spSet.nextOnBit(0); isp >= 0; isp = nprop.spSet.nextOnBit(isp + 1)) {
                final double nodeHeight = spTree.getNodeHeight(nodeID);
                {
                    double[] cp = args.cps[isp];
                    final int upi = singleStartPoints[isp];

                    int i = args.iSingle[isp];

                    while (i < cp.length && cp[i] < cHeight) {
                        if (args.indicators[upi + i] > 0) {
                            //System.out.println("  popbit s");
                            args.iSingle[isp] = i;
                            double prev = args.findPrev(cp[i], nodeHeight);
                            double mid = (prev + cp[i]) / 2.0;
                            assert nodeHeight < mid;
                            allPoints.add(new Points(mid, args.pops[upi + i]));
                        }
                        ++i;
                    }
                    args.iSingle[isp] = i;
                }

                final int kx = (isp * (2 * nSpSeqs - isp - 3)) / 2 - 1;
                for (int y = nprop.spSet.nextOnBit(isp + 1); y >= 0; y = nprop.spSet.nextOnBit(y + 1)) {

                    assert isp < y;
                    int k = kx + y;

                    double[] cp = args.cpp[k];
                    int i = args.iPair[k];
                    final int upi = pairStartPoints[k];

                    while (i < cp.length && cp[i] < cHeight) {
                        if (args.indicators[upi + i] > 0) {
                            //System.out.println("  popbit p");
                            args.iPair[k] = i;
                            final double prev = args.findPrev(cp[i], nodeHeight);
                            double mid = (prev + cp[i]) / 2.0;
                            assert nodeHeight < mid;
                            allPoints.add(new Points(mid, args.pops[upi + i]));
                        }
                        ++i;
                    }
                    args.iPair[k] = i;
                }
            }
        }

        Points[] all = null;

        if (allPoints.size() > 0) {
            all = allPoints.toArray(new Points[allPoints.size()]);
            if (all.length > 1) {
                HeapSort.sort(all);
            }

            int len = all.length;

            if (bmp) {
                int k = 0;
                while (k + 1 < len) {
                    final double t = all[k].time;
                    if (t == all[k + 1].time) {
                        int j = k + 2;
                        boolean use = all[k].use || all[k + 1].use;

                        while (j < len && t == all[j].time) {
                            use = use || all[j].use;
                            j += 1;
                        }
                        int removed = (j - k - 1);
                        all[k] = new Points(t, use);
                        for (int i = k + 1; i < len - removed; ++i) {
                            all[i] = all[i + removed];
                        }
                        len -= removed;
                    }
                    ++k;
                }
            } else {
                // duplications

                int k = 0;
                while (k + 1 < len) {
                    double t = all[k].time;
                    if (t == all[k + 1].time) {
                        int j = k + 2;
                        double v = all[k].population + all[k + 1].population;
                        while (j < len && t == all[j].time) {
                            v += all[j].population;
                            j += 1;
                        }
                        int removed = (j - k - 1);
                        all[k] = new Points(t, v / (removed + 1));
                        for (int i = k + 1; i < len - removed; ++i) {
                            all[i] = all[i + removed];
                        }
                        //System.arraycopy(all, j, all, k + 1, all.length - j + 1);
                        len -= removed;
                    }
                    ++k;
                }
            }

            if (len != all.length) {
                Points[] a = new Points[len];
                System.arraycopy(all, 0, a, 0, len);
                all = a;
            }

            if (bmp) {
                for (Points p : all) {
                    double t = p.time;
                    assert p.population == 0;
                    for (int isp = nprop.spSet.nextOnBit(0); isp >= 0; isp = nprop.spSet.nextOnBit(isp + 1)) {
                        SimpleDemographicFunction d = args.dms[isp];
                        if (t <= d.upperBound()) {
                            p.population += d.population(t);
                        }
                    }
                }
            }
        }

        pointsList[nodeID.getNumber()] = all;

        return nprop;
    }

    private int setDemographics(NodeRef nodeID, int pStart, int side, double[] pops, Points[][] pointsList) {

        final int nSpSeqs = mulspb.nSpSeqs();
        final NodeProperties nprop = props.get(nodeID);
        int pEnd;
        double p0;

        if (isExternal(nodeID)) {
            final int sps = nprop.speciesIndex;
            p0 = pops[sps];
            pEnd = pStart;
        } else {
            assert getChildCount(nodeID) == 2;

            final int iHere = setDemographics(getChild(nodeID, 0), pStart, 0, pops, pointsList);
            pEnd = setDemographics(getChild(nodeID, 1), iHere + 1, 1, pops, pointsList);
            if (constantPopulation) {
                final int i = nSpSeqs + iHere;
                p0 = pops[i];
            } else {
                final int i = nSpSeqs + iHere * 2;
                p0 = pops[i] + pops[i + 1];
            }
        }

        if (constantPopulation) {
            double[] xs = {};
            double[] ys = {p0};

            nprop.demogf = new VDdemographicFunction(xs, ys, getUnits());
            // new ConstantPopulation(p0, getUnits());
        } else {
            final double t0 = getNodeHeight(nodeID);

            Points[] p = pointsList != null ? pointsList[nodeID.getNumber()] : null;
            final int plen = p == null ? 0 : p.length;

            final boolean isRoot = nodeID == getRoot();
//        double[] xs = new double[plen + (isRoot ? 1 : 1)];
//        double[] ys = new double[plen + (isRoot ? 2 : 2)];

            final boolean useBMP = bmp && pointsList != null;
            // internal nodes add one population point for the branch end.
            // on the root (with bmp) there is no such point.
            final int len = plen + (useBMP ? (!isRoot ? 1 : 0) : 1);
            double[] xs = new double[len];
            double[] ys = new double[len + 1];
            boolean[] use = new boolean[len];
            ys[0] = p0;
            for (int i = 0; i < plen; ++i) {
                xs[i] = p[i].time - t0;
                ys[i + 1] = p[i].population;
                use[i] = p[i].use;
            }

            if (!isRoot) {
                final int anccIndex = (side == 0) ? pEnd : pStart - 1;
                final double pe = pops[nSpSeqs + anccIndex * 2 + side];
                final double b = getBranchLength(nodeID);

                xs[xs.length - 1] = b;
                ys[ys.length - 1] = pe;
            }

            if (useBMP) {
                nprop.demogf = bestLinearFit(xs, ys, use);
            } else {

                if (isRoot) {
                    // extend the last point to most ancient coalescent point. Has no effect on the demographic
                    // per se but for use when analyzing the results.

                    double h = -1;
                    for (MulSpeciesBindings.GeneTreeInfo t : mulspb.getGeneTrees()) {
                        h = Math.max(h, t.tree.getNodeHeight(t.tree.getRoot()));
                    }

                    final double rh = h - t0;
                    xs[xs.length - 1] = rh; //getNodeHeight(nodeID);
                    //spTree.setBranchLength(nodeID, rh);

                    // last value is for root branch end point
                    ys[ys.length - 1] = pointsList != null ? ys[ys.length - 2] :
                            (nonConstRootPopulation ? pops[pops.length - 1] : ys[ys.length - 2]);
                }

                nprop.demogf = new VDdemographicFunction(xs, ys, getUnits());
            }
        }
        return pEnd;
    }

    private void setNodeProperties() {
        Points[][] perBranchPoints = null;

        if (coalPointsPops != null) {
            final Args args = new Args(bmp);
            perBranchPoints = new Points[getNodeCount()][];
            getDemographicPoints(getRoot(), args, perBranchPoints);
        } else {
            // sets species info
            getDemographicPoints(getRoot(), null, null);
        }

        setDemographics(getRoot(), 0, -1, ((Parameter.Default) sppSplitPopulations).inspectParameterValues(), perBranchPoints);
    }

    private Map<NodeRef, NodeProperties> getProps() {
        if (!nodePropsReady) {
            setNodeProperties();
            nodePropsReady = true;
        }
        return props;
    }

    public DemographicFunction getNodeDemographic(NodeRef node) {
        return getProps().get(node).demogf;
    }

    public FixedBitSet spSet(NodeRef node) {
        return getProps().get(node).spSet;
    }

    public int speciesIndex(NodeRef tip) {
        assert isExternal(tip);

        // always ready even if props is dirty
        return props.get(tip).speciesIndex;
    }

    // grjtodo-oneday: not used because user starting tree not yet supported
    private Double setInitialSplitPopulations(FlexibleTree startTree, NodeRef node, int pos[]) {
        if (!startTree.isExternal(node)) {
            int loc = -1;
            for (int nc = 0; nc < startTree.getChildCount(node); ++nc) {
                final Double p = setInitialSplitPopulations(startTree, startTree.getChild(node, nc), pos);
                if (!constantPopulation && nc == 0) {
                    loc = pos[0];
                    pos[0] += 1;
                }
                if (p != null) {
                    if (constantPopulation) {
                        //
                    } else {
                        sppSplitPopulations.setParameterValueQuietly(mulspb.nSpSeqs() + 2 * loc + nc, p);
                    }
                }
            }
        }

        final String comment = (String) startTree.getNodeAttribute(node, NewickImporter.COMMENT);
        Double p0 = null;
        if (comment != null) {
            StringTokenizer st = new StringTokenizer(comment);

            p0 = Double.parseDouble(st.nextToken());
            if (startTree.isExternal(node)) {
                int ns = (Integer) startTree.getNodeAttribute(node, spIndexAttrName);
                sppSplitPopulations.setParameterValueQuietly(ns, p0);
            } else if (constantPopulation) {
                // not tested code !!
                sppSplitPopulations.setParameterValueQuietly(mulspb.nSpSeqs() + pos[0], p0);
                pos[0] += 1;
            }

            // if just one value const
            if (st.hasMoreTokens()) {
                p0 = Double.parseDouble(st.nextToken());
            }
        }
        return !constantPopulation ? p0 : null;
    }

    private SimpleTree compatibleUninformedSpeciesTree(Tree startTree) {
        double rootHeight = Double.MAX_VALUE;

        for (MulSpeciesBindings.GeneTreeInfo t : mulspb.getGeneTrees()) {
            rootHeight = Math.min(rootHeight, t.getCoalInfo()[0].ctime);
        }

        if (startTree != null) {/*
            // Allow start tree to be very basic basic - may be only partially resolved and no
            // branch lengths

            if (startTree.getExternalNodeCount() != spp.length) {
                throw new Error("Start tree error - different number of tips");
            }

            final FlexibleTree tree = new FlexibleTree(startTree, true);
            tree.resolveTree();
            final double treeHeight = tree.getRootHeight();
            if (treeHeight <= 0) {
                tree.setRootHeight(1.0);
                Utils.correctHeightsForTips(tree);
                SimpleTree.Utils.scaleNodeHeights(tree, rootHeight / tree.getRootHeight());
            }

            SimpleTree sTree = new SimpleTree(tree);
            for (int ns = 0; ns < spp.length; ns++) {
                MulSpeciesBindings.SPinfo sp = spp[ns];
                final int i = sTree.getTaxonIndex(sp.name);
                if (i < 0) {
                    throw new Error(sp.name + " is not present in the start tree");
                }
                final SimpleNode node = sTree.getExternalNode(i);
                node.setAttribute(spIndexAttrName, ns);

                // set for possible pops
                tree.setNodeAttribute(tree.getNode(tree.getTaxonIndex(sp.name)), spIndexAttrName, ns);
            }

            if (treeHeight > 0) {
                sTree.setAttribute("check", new Double(rootHeight));
            }

            {
                //assert ! constantPopulation; // not implemented yet
                int[] pos = {0};
                setInitialSplitPopulations(tree, tree.getRoot(), pos);
            }

            return sTree;*/
        }

        int nSpSeqs = mulspb.nSpSeqs();
        final double delta = rootHeight / (nSpSeqs + 1);
        double cTime = delta;

        List<SimpleNode> subs = new ArrayList<SimpleNode>(nSpSeqs);

        for (int ns = 0; ns < nSpSeqs; ns++) {
            //Taxon taxon = mulspb.getSpSeq(ns);
            final SimpleNode node = new SimpleNode();
            int sp = mulspb.spseqindex2sp(ns);
            int seq = mulspb.spseqindex2seq(ns);
            String spSeqName = mulspb.apspeciesName(sp) + seq;
            node.setTaxon(new Taxon(spSeqName));
            subs.add(node);

            node.setAttribute(spIndexAttrName, ns);
        }

        while (subs.size() > 1) {
            final SimpleNode node = new SimpleNode();
            int i = 0, j = 1;
            node.addChild(subs.get(i));
            node.addChild(subs.get(j));
            node.setHeight(cTime);
            cTime += delta;
            subs.set(j, node);
            subs.remove(i);
        }

        return new SimpleTree(subs.get(0));
    }

    public void setPreorderIndices(int[] indices) {
        setPreorderIndices(getRoot(), 0, indices);
    }

    private int setPreorderIndices(NodeRef node, int loc, int[] indices) {
        if (!isExternal(node)) {
            int l = setPreorderIndices(getChild(node, 0), loc, indices);
            indices[node.getNumber()] = l;
            loc = setPreorderIndices(getChild(node, 1), l + 1, indices);
        }
        return loc;
    }

    public String getName() {
        return getModelName();
    }

    // grj need to use a different operator
    static private MulTreeNodeSlide internalTreeOP = null;

    public int scale(double scaleFactor, int nDims, boolean testBounds) {
        assert scaleFactor > 0;
        if (nDims <= 0) {
            // actually when in an up down with operators on the gene trees the flags
            // may indicate a change
            //storeState();  // just checks assert really
            beginTreeEdit();
            final int count = getInternalNodeCount();
            for (int i = 0; i < count; ++i) {
                final NodeRef n = getInternalNode(i);
                setNodeHeight(n, getNodeHeight(n) * scaleFactor);
            }
            endTreeEdit();
            fireModelChanged(this, 1);
            return count;
        } else {
            if (nDims != 1) {
                throw new UnsupportedOperationException("not implemented for count != 1");
            }
            if (internalTreeOP == null) {
                internalTreeOP = new MulTreeNodeSlide(this, mulspb, 1);
            }

            internalTreeOP.operateOneNode(scaleFactor);
            fireModelChanged(this, 1);
            return nDims;
        }
    }

    @Override
    public boolean testBounds() {
        return true;
    }

    private final boolean verbose = false;

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (verbose) System.out.println(" SPtree: model changed " + model.getId());

        nodePropsReady = false;
        anyChange = true;
        // this should happen by default, no?
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (verbose) System.out.println(" SPtree: parameter changed " + variable.getId());

        nodePropsReady = false;
        anyChange = true;
    }

    protected void storeState() {
        assert !treeChanged;
        assert !anyChange;
    }

    protected void restoreState() {
        if (verbose) System.out.println(" SPtree: restore (" + treeChanged + "," + anyChange + ")");

        if (treeChanged) {
            //
            spTree.beginTreeEdit();

            for (int k = 0; k < getInternalNodeCount(); ++k) {
                final NodeRef node = getInternalNode(k);
                final int index = node.getNumber();
                final double h = heights[index];
                if (getNodeHeight(node) != h) {
                    setNodeHeight(node, h);
                }
                for (int nc = 0; nc < 2; ++nc) {
                    final NodeRef child = getChild(node, nc);

                    final NodeRef child1 = children[2 * index + nc];
                    if (child != child1) {
                        replaceChild(node, child, child1);
                    }
                    assert getParent(child1) == node;
                }
            }
            setRoot(children[children.length - 1]);

            if (verbose) System.out.println("  restored to: " + spTree);

            spTree.endTreeEdit();
        }
        if (treeChanged || anyChange) {
            setNodeProperties();
        }
        treeChanged = false;
        anyChange = false;
    }

    protected void acceptState() {
        if (verbose) System.out.println(" SPtree: accept");

        treeChanged = false;
        anyChange = false;
    }

    String previousTopology = null;

    public boolean logNow(long state) {
        final String curTop = TreeUtils.uniqueNewick(spTree, spTree.getRoot());
        if (state == 0 || !curTop.equals(previousTopology)) {
            previousTopology = curTop;
            return true;
        }
        return false;
    }

//    TreeTrait dmf = new TreeTrait.S() {
//        public String getTraitName() {
//            return "dmf";
//        }
//
//        public Intent getIntent() {
//            return Intent.NODE;
//        }
//
//        public String getTrait(Tree tree, NodeRef node) {
//            assert tree == SpeciesTreeModel.this;
//
//            //final VDdemographicFunction df = getProps().get(node).demogf;
//
//            final DemographicFunction df = getNodeDemographic(node);
//            return df.toString();
//        }
//    };

    TreeTrait dmt = new TreeTrait.DA() {
        public String getTraitName() {
            return "dmt";
        }

        public Intent getIntent() {
            return Intent.NODE;
        }

        public double[] getTrait(Tree tree, NodeRef node) {
            assert tree == MulSpeciesTreeModel.this;

            final VDdemographicFunction df = (VDdemographicFunction) getNodeDemographic(node);
            return df.times();
        }
    };

    TreeTrait dmv = new TreeTrait.DA() {
        public String getTraitName() {
            return "dmv";
        }

        public Intent getIntent() {
            return Intent.NODE;
        }

        public double[] getTrait(Tree tree, NodeRef node) {
            assert tree == MulSpeciesTreeModel.this;

            final VDdemographicFunction df = (VDdemographicFunction) getNodeDemographic(node);
            return df.values();
        }
    };


    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[]{dmt, dmv};
    }

    public TreeTrait getTreeTrait(String key) {
        if (key.equals(dmt.getTraitName())) {
            return dmt;
        } else if (key.equals(dmv.getTraitName())) {
            return dmv;
        }
        throw new IllegalArgumentException();
    }


    // boring delegation

    public SimpleTree getSimpleTree() {
        return spTree;
    }

    public Tree getCopy() {
        return spTree.getCopy();
    }

    public Type getUnits() {
        return spTree.getUnits();
    }

    public void setUnits(Type units) {
        spTree.setUnits(units);
    }

    public int getNodeCount() {
        return spTree.getNodeCount();
    }

    public boolean hasNodeHeights() {
        return spTree.hasNodeHeights();
    }

    public double getNodeHeight(NodeRef node) {
        return spTree.getNodeHeight(node);
    }

    public double getNodeRate(NodeRef node) {
        return spTree.getNodeRate(node);
    }

    public Taxon getNodeTaxon(NodeRef node) {
        return spTree.getNodeTaxon(node);
    }

    public int getChildCount(NodeRef node) {
        return spTree.getChildCount(node);
    }

    public boolean isExternal(NodeRef node) {
        return spTree.isExternal(node);
    }

    public boolean isRoot(NodeRef node) {
        return spTree.isRoot(node);
    }

    public NodeRef getChild(NodeRef node, int i) {
        return spTree.getChild(node, i);
    }

    public NodeRef getParent(NodeRef node) {
        return spTree.getParent(node);
    }

    public boolean hasBranchLengths() {
        return spTree.hasBranchLengths();
    }

    public double getBranchLength(NodeRef node) {
        return spTree.getBranchLength(node);
    }

    public void setBranchLength(NodeRef node, double length) {
        spTree.setBranchLength(node, length);
    }

    public NodeRef getExternalNode(int i) {
        return spTree.getExternalNode(i);
    }

    public NodeRef getInternalNode(int i) {
        return spTree.getInternalNode(i);
    }

    public NodeRef getNode(int i) {
        return spTree.getNode(i);
    }

    public int getExternalNodeCount() {
        return spTree.getExternalNodeCount();
    }

    public int getInternalNodeCount() {
        return spTree.getInternalNodeCount();
    }

    public NodeRef getRoot() {
        return spTree.getRoot();
    }

    public void setRoot(NodeRef r) {
        spTree.setRoot(r);
    }

    public void addChild(NodeRef p, NodeRef c) {
        spTree.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {
        spTree.removeChild(p, c);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        spTree.replaceChild(node, child, newChild);
    }

    public boolean beginTreeEdit() {
        boolean beingEdited = spTree.beginTreeEdit();
        if (!beingEdited) {
            // save tree for restore
            for (int n = 0; n < getInternalNodeCount(); ++n) {
                final NodeRef node = getInternalNode(n);
                final int k = node.getNumber();
                children[2 * k] = getChild(node, 0);
                children[2 * k + 1] = getChild(node, 1);
                heights[k] = getNodeHeight(node);
            }
            children[children.length - 1] = getRoot();

            treeChanged = true;
            nodePropsReady = false;
            //anyChange = true;
        }
        return beingEdited;
    }

    public void endTreeEdit() {
        spTree.endTreeEdit();
        fireModelChanged();
    }

    public void setNodeHeight(NodeRef n, double height) {
        spTree.setNodeHeight(n, height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        spTree.setNodeRate(n, rate);
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        spTree.setNodeAttribute(node, name, value);
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        return spTree.getNodeAttribute(node, name);
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        return spTree.getNodeAttributeNames(node);
    }

    public int getTaxonCount() {
        return spTree.getTaxonCount();
    }

    public Taxon getTaxon(int taxonIndex) {
        return spTree.getTaxon(taxonIndex);
    }

    public String getTaxonId(int taxonIndex) {
        return spTree.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(String id) {
        return spTree.getTaxonIndex(id);
    }

    public int getTaxonIndex(Taxon taxon) {
        // can't compare taxa
        return getTaxonIndex(taxon.getId());
    }

    public List<Taxon> asList() {
        return spTree.asList();
    }

    public Iterator<Taxon> iterator() {
        return spTree.iterator();
    }

    public Object getTaxonAttribute(int taxonIndex, String name) {
        return spTree.getTaxonAttribute(taxonIndex, name);
    }

    public int addTaxon(Taxon taxon) {
        return spTree.addTaxon(taxon);
    }

    public boolean removeTaxon(Taxon taxon) {
        return spTree.removeTaxon(taxon);
    }

    public void setTaxonId(int taxonIndex, String id) {
        spTree.setTaxonId(taxonIndex, id);
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        spTree.setTaxonAttribute(taxonIndex, name, value);
    }

    public String getId() {
        return spTree.getId();
    }

    public void setId(String id) {
        spTree.setId(id);
    }

    public void setAttribute(String name, Object value) {
        spTree.setAttribute(name, value);
    }

    public Object getAttribute(String name) {
        return spTree.getAttribute(name);
    }

    public Iterator<String> getAttributeNames() {
        return spTree.getAttributeNames();
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        spTree.addMutableTreeListener(listener);
    }

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        spTree.addMutableTaxonListListener(listener);
    }

    public static Parameter createCoalPointsPopParameter(MulSpeciesBindings spb, Double value, Boolean bmp) {
        int dim = 0;
        for (double[] d : spb.getPopTimesSingle()) {
            dim += d.length;
        }

        if (!bmp) {
            for (double[] d : spb.getPopTimesPair()) {
                dim += d.length;
            }
        }

        return new Parameter.Default(dim, value);
    }

    public static Parameter createSplitPopulationsParameter(MulSpeciesBindings spb, double value, boolean root, boolean constPop) {
        int dim;
        if (constPop) {
            // one per node
            dim = 2 * spb.nSpSeqs() - 1;
        } else {
            // one per species leaf (ns) + 2 per internal node (2*(ns-1)) + optionally one for the root
            dim = 3 * spb.nSpSeqs() - 2 + (root ? 1 : 0);
        }
        return new Parameter.Default(dim, value);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SPECIES_MODELS;
    }

    @Override
    public String getDescription() {
        return "Multiply labelled species tree";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(
                new Citation(
                        new Author[]{
                                new Author("GR", "Jones")
                        },
                        Citation.Status.IN_PREPARATION
                ));
    }

}
