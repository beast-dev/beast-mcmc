/*
 * PopsIOSpeciesTreeModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.speciation;

import java.util.*;
import java.util.logging.Logger;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.evomodelxml.speciation.PopsIOSpeciesTreeModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.Scalable;
import dr.util.AlloppMisc;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import jebl.util.FixedBitSet;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

/**
 * @author  Graham  Jones
 * Date: 10/05/2012
 */

/*
nodes[] and root are the fundamental data here. They are used for
the SlidableTree implementation, not just supporting the MCMC move,
but also the calculations dealing with coalescences (compatibility
and likelihood).

oldtree stores the state so it can be restored after a rejected move. It
only stores taxa at tips, topology and node times. The other data (unions, lineages,
coalescences) are reconstructed from that.

stree is used to implement Tree, which is used for logging and for
calculation of the tree prior (eg Yule model). stree is a copy of
topology and node times  from nodes[].
*/


public class PopsIOSpeciesTreeModel extends AbstractModel implements SlidableTree, Tree,
                                               Scalable, TreeLogger.LogUpon, Citable {
    private final PopsIOSpeciesBindings piosb;
    private final Parameter popPriorScale;
    private final PriorComponent[] priorComponents;
    private PopsIONode[] pionodes;
    private int rootn;
    private PopsIONode[] oldpionodes;
    private int oldrootn;
    private SimpleTree stree;

    @Override
    public boolean logNow(long state) {
        // for debugging, set logEvery=0 in XML:
        //    <!--  species tree log file.  -->
        //    <logTree id="pioTreeFileLog" logEvery="0" fileName="C:\Users\....

        if (state % 10000 == 0) {
            System.out.println("DEBUGGING: PopsIOSpeciesTreeModel.logNow(), state = " + state);
        }
        return (state % 2500) == 0;
    }


    public static class PriorComponent {
        private final double weight;
        private final double alpha;
        private final double beta;

        // inv gamma pdf is parameterized as  b^a/Gamma(a)  x^(-a-1)  exp(-b/x)
        // mean is b/(a-1) if a>1, var is  b^2/((a-1)^2 (a-2)) if a>2.

        public PriorComponent(double weight, double alpha, double beta) {
            this.weight = weight;
            this.alpha = alpha;
            this.beta = beta;
        }
    }


    /*
     *  The parent, child, height, taxon fileds implement the basic binary tree.
     *  The other fields are 'working' fields.
     *
     *  nodeNumber is required to implement NodeRef, needed for replaceSlidableRoot()
     *  (possibly one can manage this another way?)
     *
     *  union, coalheights, nlineages are used for calculations: determining compatibility
     *  with gene trees and calculating likelihood.
     */
    private class PopsIONode extends AlloppNode.Abstract implements AlloppNode, NodeRef {
        private int anc;
        private int lft;
        private int rgt;
        private double height;
        private Taxon taxon;
        private FixedBitSet union;
        private ArrayList<Double> coalheights;
        private int nlineages;
        private int coalcount;
        private double coalintensity;
        private int nodeNumber;

        // dud constuctor
        PopsIONode(int nn) {
            anc = -1;
            lft = -1;
            rgt = -1;
            height = -1.0;
            coalheights = new ArrayList<Double>();
            taxon = new Taxon("");
            union = null;
            coalcount = 0;
            coalintensity = 0.0;
            nodeNumber = nn;
        }



        // copy constructor
        public PopsIONode(PopsIONode node) {
            anc = node.anc;
            lft = node.lft;
            rgt = node.rgt;
            nodeNumber = node.nodeNumber;
            copyNonTopologyFields(node);
        }




        private void copyNonTopologyFields(PopsIONode node) {
            height = node.height;
            taxon = new Taxon(node.taxon.getId());
            nlineages = node.nlineages;
            if (node.union == null) {
                union = null;
            } else {
                union = new FixedBitSet(node.union);
            }
            coalheights = new ArrayList<Double>();
            for (int i = 0; i < node.coalheights.size(); i++) {
                coalheights.add(node.coalheights.get(i));
            }
        }




        public String asText(int indentlen) {
            StringBuilder s = new StringBuilder();
            Formatter formatter = new Formatter(s, Locale.US);
            if (lft < 0) {
                formatter.format("%s ", taxon.getId());
            } else {
                formatter.format("%s ", "+");
            }
            while (s.length() < 20-indentlen) {
                formatter.format("%s", " ");
            }
            formatter.format("%s ", AlloppMisc.nonnegIn8Chars(height));
            formatter.format("%20s ", AlloppMisc.FixedBitSetasText(union));
            formatter.format("%3d  ", nlineages);
            for (int c = 0; c < coalheights.size(); c++) {
                formatter.format(AlloppMisc.nonnegIn8Chars(coalheights.get(c)) + ",");
            }
            return s.toString();
        }


        @Override
        public int getNumber() {
            return nodeNumber;
        }

        @Override
        public void setNumber(int n) {
            nodeNumber = n;
        }

        @Override
        public int nofChildren() {
            return (lft < 0) ? 0 : 2;
        }

        @Override
        public AlloppNode getChild(int ch) {
            return ch==0 ? pionodes[lft] : pionodes[rgt];
        }

        @Override
        public AlloppNode getAnc() {
            return pionodes[anc];
        }

        @Override
        public Taxon getTaxon() {
            return taxon;
        }

        @Override
        public double getHeight() {
            return height;
        }

        @Override
        public FixedBitSet getUnion() {
            return union;
        }

        @Override
        public void setChild(int ch, AlloppNode newchild) {
            int newch = ((PopsIONode)newchild).nodeNumber;
            if (ch == 0) {
                lft = newch;
            } else {
                rgt = newch;
            }
        }

        @Override
        public void setAnc(AlloppNode anc) {
            this.anc = ((PopsIONode)anc).nodeNumber;

        }

        @Override
        public void setTaxon(String name) {
            this.taxon = new Taxon(name);

        }

        @Override
        public void setHeight(double height) {
            this.height = height;

        }

        @Override
        public void setUnion(FixedBitSet union) {
            this.union = union;

        }

        @Override
        public void addChildren(AlloppNode c0, AlloppNode c1) {
            lft = ((PopsIONode)c0).nodeNumber;
            pionodes[lft].anc = nodeNumber;
            rgt = ((PopsIONode)c1).nodeNumber;
            pionodes[rgt].anc = nodeNumber;
        }

    }


    public PopsIOSpeciesTreeModel(PopsIOSpeciesBindings piosb,  Parameter popPriorScale, PriorComponent[] priorComponents) {
        super(PopsIOSpeciesTreeModelParser.PIO_SPECIES_TREE);
        this.piosb = piosb;

        this.popPriorScale = popPriorScale;
        addVariable(popPriorScale);
        popPriorScale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.priorComponents = priorComponents;

        PopsIOSpeciesBindings.SpInfo[] species = piosb.getSpecies();
        int nTaxa = species.length;
        int nNodes = 2 * nTaxa - 1;
        pionodes = new PopsIONode[nNodes];
        for (int n = 0; n < nNodes; n++) {
            pionodes[n] = new PopsIONode(n);
        }
        ArrayList<Integer> tojoin = new ArrayList<Integer>(nTaxa);
        for (int n = 0; n < nTaxa; n++) {
            pionodes[n].setTaxon(species[n].name);
            pionodes[n].setHeight(0.0);
            pionodes[n].setUnion(piosb.tipUnionFromTaxon(pionodes[n].getTaxon()));
            tojoin.add(n);
        }
        double rate = 1.0;
        double treeheight = 0.0;
        for (int i = 0; i < nTaxa-1; i++) {
            int numtojoin = tojoin.size();
            int j = MathUtils.nextInt(numtojoin);
            Integer child0 = tojoin.get(j);
            tojoin.remove(j);
            int k = MathUtils.nextInt(numtojoin-1);
            Integer child1 = tojoin.get(k);
            tojoin.remove(k);
            pionodes[nTaxa+i].addChildren(pionodes[child0],pionodes[child1]);
            pionodes[nTaxa+i].setHeight(treeheight + randomnodeheight(numtojoin*rate));
            treeheight = pionodes[nTaxa+i].getHeight();
            tojoin.add(nTaxa+i);
        }
        rootn = pionodes.length - 1;

        double scale = 0.99 * piosb.initialMinGeneNodeHeight() / pionodes[rootn].height;
        scaleAllHeights(scale);
        pionodes[rootn].fillinUnionsInSubtree(piosb.getSpecies().length);

        stree = makeSimpleTree();

        Logger.getLogger("dr.evomodel.speciation.popsio").info("\tConstructing a PopsIO Species Tree Model, please cite:\n"
                + Citable.Utils.getCitationString(this));
    }



    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(
                new Author[]{
                        new Author("Graham", "Jones")
                },
                "WORKING TITLE: A multi-species coalescent model with population parameters integrated out",
                "??",  // journal
                Citation.Status.IN_PREPARATION
        ));
        return citations;
    }




    public String toString() {
        int ngt = piosb.numberOfGeneTrees();
        String nl = System.getProperty("line.separator");
        String s = nl + pioTreeAsText() + nl;
        for (int g = 0; g < ngt; g++) {
            s += "Gene tree " + g + nl;
            s += piosb.genetreeAsText(g) + nl;
        }
        s += nl;
        return s;
    }



    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[1];
        columns[0] = new LogColumn.Default("    species-tree and gene trees", this);
        return columns;
    }



    private int scaleAllHeights(double scale) {
        for (int nn = 0; nn < pionodes.length; nn++) {
            pionodes[nn].height *= scale;
        }
        return pionodes.length;
    }


    String pioTreeAsText() {
        String header = "topology             height         union         nlin coalheights" + System.getProperty("line.separator");

        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + subtreeAsText(pionodes[rootn], s, x, 0, "");
    }

    /*
      * Called from PopsIOSpeciesBindings to check if a node in a gene tree
      * is compatible with the network.
      */
    public boolean coalescenceIsCompatible(double height, FixedBitSet union) {
        PopsIONode node = (PopsIONode) pionodes[rootn].nodeOfUnionInSubtree(union);
        return (node.height <= height);
    }



    /*
    * Called from PopsIOSpeciesBindings to zero coalescent information
    * in branches of species tree. Called once for all genes.
    */
    public void zeroCoalCountsIntensities() {
        zeroSubtreeCoalCountsIntensities(pionodes[rootn]);
    }

    /*
    * Called from PopsIOSpeciesBindings to remove gene-specific coalescent information
    * from branches of species tree. Required before call to recordCoalescence.
    */
    public void clearCoalescences() {
        clearSubtreeCoalescences(pionodes[rootn]);
    }



    /*
      * Called from PopsIOSpeciesBindings to add a node from a gene tree
      * to its branch in species tree.
      */
    public void recordCoalescence(double height, FixedBitSet union) {
        PopsIONode node = (PopsIONode) pionodes[rootn].nodeOfUnionInSubtree(union);
        assert (node.height <= height);
        while (node.anc >= 0  &&  pionodes[node.anc].height <= height) {
            node = pionodes[node.anc];
        }
        node.coalheights.add(height);

    }

    /*
    * Called from PopsIOSpeciesBindings to sort coalecent times of a gene
    * in the branch of a node in species tree.
    */
    public void sortCoalescences() {
        for (PopsIONode node : pionodes) {
            Collections.sort(node.coalheights);
        }
    }

    /*
      * Called from PopsIOSpeciesBindings to record the number
      * of gene lineages at nodes of species tree.
      */
    public void recordLineageCounts() {
        recordSubtreeLineageCounts(pionodes[rootn]);
    }


    /*
      * Called from PopsIOSpeciesBindings to accumulate information
      * for all genes at nodes of species tree.
      */
    public void accumCoalCountsIntensities() {
        accumSubtreeCoalCountsIntensities(pionodes[rootn]);
    }


    public double logLhoodAllGeneTreesInSpeciesTree() {
        double [] llhoodcpts = new double[priorComponents.length];

        double totalweight = 0.0;
        for (int i = 0; i < priorComponents.length; i++) {
            totalweight += priorComponents[i].weight;
        }
        for (int i = 0; i < priorComponents.length; i++) {
            llhoodcpts[i] = Math.log(priorComponents[i].weight / totalweight);
            double sigma = popPriorScale.getParameterValue(0);
            double alpha = priorComponents[i].alpha;
            double beta = sigma * priorComponents[i].beta;
            llhoodcpts[i] += logLhoodAllGeneTreesInSpeciesSubtree(pionodes[rootn], alpha, beta);
        }
        return logsumexp(llhoodcpts);
    }



    private double logsumexp(double x[]) {
        double maxx = Double.MIN_VALUE;
        for (double d : x) {
            if (d > maxx) { maxx = d; }
        }
        double sum = 0.0;
        for (double d : x) {
            sum += Math.exp(d-maxx);
        }
        return maxx + Math.log(sum);
    }


    public void fixupAfterNodeSlide() {
        pionodes[rootn].fillinUnionsInSubtree(piosb.getSpecies().length);
        stree = makeSimpleTree();
    }




    private String subtreeAsText(PopsIONode node, String s, Stack<Integer> x, int depth, String b) {
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
        s += node.asText(indent.length());
        s += System.getProperty("line.separator");
        String subs = "";
        if (node.lft >= 0) {
            x.push(depth);
            subs += subtreeAsText(pionodes[node.lft], "", x, depth+1, "-");
            x.pop();
            subs += subtreeAsText(pionodes[node.rgt], "", x, depth+1, "`-");
        }
        return s + subs;
    }


    private double logLhoodAllGeneTreesInSpeciesSubtree(PopsIONode node, double alpha, double beta) {
        double loglike = 0.0;
        if (node.lft >= 0) {
            loglike += logLhoodAllGeneTreesInSpeciesSubtree(pionodes[node.lft], alpha, beta);
            loglike += logLhoodAllGeneTreesInSpeciesSubtree(pionodes[node.rgt], alpha, beta);
        }
        loglike += branchLLInPopsIOtree(node, alpha, beta);
        return loglike;
    }


    private double branchLLInPopsIOtree(PopsIONode node, double alpha, double beta) {
        int q = node.coalcount;
        double gamma = node.coalintensity;
        double llhood = 0.0;
        for (int i = 1; i <= q-1; i++) {
            llhood += Math.log(alpha+i);
        }
        llhood += alpha * Math.log(beta);
        llhood -= (alpha + q) * Math.log(beta + gamma);
        return llhood;
    }

    private SimpleTree makeSimpleTree() {
        SimpleNode[] snodes = new SimpleNode[pionodes.length];
        for (int n = 0; n < pionodes.length; n++) {
            snodes[n] = new SimpleNode();
            snodes[n].setTaxon(null);
        }
        makesimplesubtree(snodes, 0, pionodes[rootn]);
        return new SimpleTree(snodes[pionodes.length-1]);
    }


    // for makeSimpleTree()
    private int makesimplesubtree(SimpleNode[] snodes, int nextsn, PopsIONode pionode) {
        if (pionode.lft < 0) {
            Taxon tx = new Taxon(pionode.taxon.getId());
            if (nextsn >= snodes.length) {
                System.out.println("BUG: makesimplesubtree()");
            }
            snodes[nextsn].setTaxon(tx);
        } else {
            nextsn = makesimplesubtree(snodes, nextsn, pionodes[pionode.lft]);
            int subtree0 = nextsn-1;
            nextsn = makesimplesubtree(snodes, nextsn, pionodes[pionode.rgt]);
            int subtree1 = nextsn-1;
            snodes[nextsn].addChild(snodes[subtree0]);
            snodes[nextsn].addChild(snodes[subtree1]);
        }
        snodes[nextsn].setHeight(pionode.height);
        return nextsn+1;
    }



    private void zeroSubtreeCoalCountsIntensities(PopsIONode node) {
        if (node.lft >= 0) {
            zeroSubtreeCoalCountsIntensities(pionodes[node.lft]);
            zeroSubtreeCoalCountsIntensities(pionodes[node.rgt]);
        }
        node.coalcount = 0;
        node.coalintensity = 0.0;
    }



    private void clearSubtreeCoalescences(PopsIONode node) {
        if (node.lft >= 0) {
            clearSubtreeCoalescences(pionodes[node.lft]);
            clearSubtreeCoalescences(pionodes[node.rgt]);
        }
        if (node == null) {
            System.out.println("BUG");
        }
        if (node.coalheights == null) {
            System.out.println("BUG");
        }
        node.coalheights.clear();
    }


    private void recordSubtreeLineageCounts(PopsIONode node) {
        if (node.lft < 0) {
            int spIndex = piosb.speciesId2index(node.getTaxon().getId());
            node.nlineages = piosb.nLineages(spIndex);
        } else {
            node.nlineages = 0;
            recordSubtreeLineageCounts(pionodes[node.lft]);
            node.nlineages += pionodes[node.lft].nlineages - pionodes[node.lft].coalheights.size();
            recordSubtreeLineageCounts(pionodes[node.rgt]);
            node.nlineages += pionodes[node.rgt].nlineages - pionodes[node.rgt].coalheights.size();
        }
    }


    private void accumSubtreeCoalCountsIntensities(PopsIONode node) {
        if (node.lft >= 0) {
            accumSubtreeCoalCountsIntensities(pionodes[node.lft]);
            accumSubtreeCoalCountsIntensities(pionodes[node.rgt]);
        }
        int k = node.coalheights.size();
        node.coalcount += k;
        double [] t = new double[k + 2];
        t[0] = node.height;
        for (int i = 0; i < k; i++) {
            t[i + 1] = node.coalheights.get(i);
        }
        t[k+1] = (node.anc < 0) ? piosb.maxGeneTreeHeight() : pionodes[node.anc].height;
        int n = node.nlineages;
        for (int i = 0; i <= k; i++) {
            node.coalintensity += (t[i + 1] - t[i]) * 0.5 * (n - i) * (n - i - 1);
        }
    }


    private double randomnodeheight(double rate) {
        return MathUtils.nextExponential(rate) + 1e-6/rate;
        // 1e-6/rate to avoid very tiny heights
    }



    /*******************************************************************************/
    // for Scalable.
    /*******************************************************************************/
    @Override
    public int scale(double factor, int nDims) throws OperatorFailedException {
        int n = scaleAllHeights(factor);
        stree = makeSimpleTree();
        return n;
    }


    @Override
    public String getName() {
        return PopsIOSpeciesTreeModelParser.PIO_SPECIES_TREE;

    }
    /*******************************************************************************/
    // for AbstractModel.
    /*******************************************************************************/
    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    @Override
    protected void storeState() {
        oldpionodes = new PopsIONode[pionodes.length];
        for (int n = 0; n < oldpionodes.length; n++) {
            oldpionodes[n] = new PopsIONode(pionodes[n]);
        }
        oldrootn = rootn;
    }

    @Override
    protected void restoreState() {
        pionodes = new PopsIONode[oldpionodes.length];
        for (int n = 0; n < pionodes.length; n++) {
            pionodes[n] = new PopsIONode(oldpionodes[n]);
        }
        rootn = oldrootn;
        stree = makeSimpleTree();
    }

    @Override
    protected void acceptState() {
    }



    /*******************************************************************************/
    // for SlidableTree.
    /*******************************************************************************/
    @Override
    public NodeRef getSlidableRoot() {
        assert pionodes[rootn].anc < 0;
        return pionodes[rootn];

    }

    @Override
    public void replaceSlidableRoot(NodeRef newroot) {
        rootn = newroot.getNumber();
        pionodes[rootn].anc = -1;
    }

    @Override
    public int getSlidableNodeCount() {
        return pionodes.length;
    }

    @Override
    public Taxon getSlidableNodeTaxon(NodeRef node) {
        assert node == pionodes[node.getNumber()];
        return ((PopsIONode)node).getTaxon();
    }


    @Override
    public double getSlidableNodeHeight(NodeRef node) {
        assert node == pionodes[node.getNumber()];
        return ((PopsIONode)node).getHeight();
    }

    @Override
    public void setSlidableNodeHeight(NodeRef node, double height) {
        assert node == pionodes[node.getNumber()];
        ((PopsIONode)node).height = height;
    }

    @Override
    public boolean isExternalSlidable(NodeRef node) {
        return (pionodes[node.getNumber()].lft < 0);
    }

    @Override
    public NodeRef getSlidableChild(NodeRef node, int j) {
        int n = node.getNumber();
        return j == 0 ? pionodes[ pionodes[n].lft ] : pionodes[ pionodes[n].rgt ];
    }

    @Override
    public void replaceSlidableChildren(NodeRef node, NodeRef lft, NodeRef rgt) {
        int nn = node.getNumber();
        int lftn = lft.getNumber();
        int rgtn = rgt.getNumber();
        assert pionodes[nn].lft >= 0;
        pionodes[nn].lft = lftn;
        pionodes[nn].rgt = rgtn;
        pionodes[lftn].anc = pionodes[nn].nodeNumber;
        pionodes[rgtn].anc = pionodes[nn].nodeNumber;
    }



    /*******************************************************************************/
    // For Tree
    /*******************************************************************************/

    @Override
    public NodeRef getRoot() {
        return stree.getRoot();
    }

    @Override
    public int getNodeCount() {
        return stree.getNodeCount();
    }

    @Override
    public NodeRef getNode(int i) {
        return stree.getNode(i);
    }

    @Override
    public NodeRef getInternalNode(int i) {
        return stree.getInternalNode(i);
    }

    @Override
    public NodeRef getExternalNode(int i) {
        return stree.getExternalNode(i);
    }

    @Override
    public int getExternalNodeCount() {
        return stree.getExternalNodeCount();
    }

    @Override
    public int getInternalNodeCount() {
        return stree.getInternalNodeCount();
    }

    @Override
    public Taxon getNodeTaxon(NodeRef node) {
        return stree.getNodeTaxon(node);
    }

    @Override
    public boolean hasNodeHeights() {
        return stree.hasNodeHeights();
    }

    @Override
    public double getNodeHeight(NodeRef node) {
        return stree.getNodeHeight(node);
    }

    @Override
    public boolean hasBranchLengths() {
        return stree.hasBranchLengths();
    }

    @Override
    public double getBranchLength(NodeRef node) {
        return stree.getBranchLength(node);
    }

    @Override
    public double getNodeRate(NodeRef node) {
        return stree.getNodeRate(node);
    }

    @Override
    public Object getNodeAttribute(NodeRef node, String name) {
        return stree.getNodeAttribute(node, name);
    }

    @Override
    public Iterator getNodeAttributeNames(NodeRef node) {
        return stree.getNodeAttributeNames(node);
    }

    @Override
    public boolean isExternal(NodeRef node) {
        return stree.isExternal(node);
    }

    @Override
    public boolean isRoot(NodeRef node) {
        return stree.isRoot(node);
    }

    @Override
    public int getChildCount(NodeRef node) {
        return stree.getChildCount(node);
    }

    @Override
    public NodeRef getChild(NodeRef node, int j) {
        return stree.getChild(node, j);
    }

    @Override
    public NodeRef getParent(NodeRef node) {
        return stree.getParent(node);
    }

    @Override
    public Tree getCopy() {
        return stree.getCopy();
    }

    @Override
    public void setAttribute(String name, Object value) {
        stree.setAttribute(name, value);
    }

    @Override
    public Object getAttribute(String name) {
        return stree.getAttribute(name);
    }

    @Override
    public Iterator<String> getAttributeNames() {
        return stree.getAttributeNames();
    }

    @Override
    public int getTaxonCount() {
        return stree.getTaxonCount();
    }

    @Override
    public Taxon getTaxon(int taxonIndex) {
        return stree.getTaxon(taxonIndex);
    }

    @Override
    public String getTaxonId(int taxonIndex) {
        return stree.getTaxonId(taxonIndex);
    }

    @Override
    public int getTaxonIndex(String id) {
        return stree.getTaxonIndex(id);
    }

    @Override
    public int getTaxonIndex(Taxon taxon) {
        return stree.getTaxonIndex(taxon);
    }

    @Override
    public List<Taxon> asList() {
        return stree.asList();
    }

    @Override
    public Object getTaxonAttribute(int taxonIndex, String name) {
        return stree.getTaxonAttribute(taxonIndex, name);
    }

    @Override
    public Iterator<Taxon> iterator() {
        return stree.iterator();
    }

    @Override
    public Type getUnits() {
        return stree.getUnits();
    }

    @Override
    public void setUnits(Type units) {
        stree.setUnits(units);
    }


}

