

/*
 * AlloppMulLabTree.java
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Formatter;
import java.util.Locale;
import java.util.Stack;

import dr.evolution.tree.*;
import jebl.util.FixedBitSet;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.evomodel.alloppnet.util.AlloppMisc;





/**
 * An AlloppMulLabTree represents the species network as single
 * binary tree with tips that can be multiply labelled with species.
 *
 * @author Graham Jones
 *         Date: 13/09/2011
 */


/*
 * An AlloppMulLabTree represents the species network as single
 * binary tree with tips that can be multiply labelled with species.
 *
 * classes LegLink and FootLinks are for gathering and organising
 * the links between trees of different ploidy, so that the
 * rootward-pointing legs can become tipward-pointing branches.
 *
 * SpSqUnion is used for sorting the nodes in an AlloppMulLabTree. It is
 * used by Comparator SPUNION_ORDER, and hence indirectly by
 * fillinpopvals().
 *
 * class BranchPopulationAndLineages records the information needed
 * to calculate the probability of coalescences in a single branch of the
 * AlloppMulLabTree.
 *
 */


/* mlnodes[], rootn implement the tree; nextn is for building it
 *
 * apsp references the (species, indivs, sequences) structure.
 *
 * popvals references the population parameters. fillinpopvals() assigns
 * them to branches.
 *
 * simptree is so that AlloppSpeciesNetworkModel, which contains a AlloppMulLabTree,
 * can implement the Tree interface.
 */
public class AlloppMulLabTree  {
    private MulLabNode[] mlnodes;
    private int rootn;
    private int nextn;
    private AlloppSpeciesBindings apsp;
    private Parameter tippopvals;
    private Parameter rootpopvals;
    private double [] hybpopvals;
    private int nofhybpopvals;
    public SimpleTree simptree;




    /*
    * parent, child[] join the nodes into a binary tree.
    *
    * height is time into past
    *
    * popsize is population at start of branch, and for tips, also at end.
    *
    * union is a set of species for a single choice of sequence copy from
    * each individual of the species. There is one bit for each Taxon like
    * "c0" or "c1"
    */
    private class MulLabNode extends AlloppNode.Abstract implements AlloppNode, NodeRef {
        private int nodeNumber;
        private int anc;
        private int lft;
        private int rgt;
        private double height;
        private boolean tetraroot;   // true iff this node is root of a teraploid subtree
        private boolean intetratree;  // true iff this is in a tetratree
        private boolean tetraancestor;  // true iff this is a node with no diploid tip descendants
        private int ttreeindex;
        private double hybridheight;
        private FixedBitSet union;
        private ArrayList<Double> coalheights;
        private int nlineages;
        private Taxon taxon;
        private int tippopindex;
        private int hybpopindex;
        private int rootpopindex;

        // dud constuctor
        MulLabNode(int nn) {
            nodeNumber = nn;
            anc = -1;
            lft = -1;
            rgt = -1;
            height = -1.0;
            tetraroot = false;
            intetratree = false;
            tetraancestor = false;
            ttreeindex = -1;
            hybridheight = -1.0;
            coalheights = new ArrayList<Double>();
            taxon = new Taxon("");
            tippopindex = -1;
            hybpopindex	= -1;
            rootpopindex = -1;
        }


        public double tippop() {
            return tippopvals.getParameterValue(tippopindex);
        }

        public double hybpop() {
            return hybpopvals[hybpopindex];
        }

        public double rootpop() {
            return rootpopvals.getParameterValue(rootpopindex);
        }

        @Override
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
            formatter.format("%60s ", AlloppMisc.FixedBitSetasText(union));
            double tippop = (tippopindex >= 0) ? tippop() : -1.0;
            formatter.format("%s %s ", AlloppMisc.nonnegIntIn2Chars(tippopindex), AlloppMisc.nonnegIn8Chars(tippop));
            double hybpop = (hybpopindex >= 0) ? hybpop() : -1.0;
            formatter.format("%s %s ", AlloppMisc.nonnegIntIn2Chars(hybpopindex), AlloppMisc.nonnegIn8Chars(hybpop));
            double rootpop = (rootpopindex >= 0) ? rootpop() : -1.0;
            formatter.format("%s %s ", AlloppMisc.nonnegIntIn2Chars(rootpopindex), AlloppMisc.nonnegIn8Chars(rootpop));
            formatter.format("%s ", tetraroot ? "tetroot" : "       ");
            formatter.format("%s ", AlloppMisc.nonnegIn8Chars(hybridheight));
            formatter.format("%3d  ", nlineages);
            for (int c = 0; c < coalheights.size(); c++) {
                formatter.format(AlloppMisc.nonnegIn8Chars(coalheights.get(c)) + ",");
            }
            return s.toString();
        }


        @Override
        public int nofChildren() {
            return ((lft < 0) ? 0 : 2);
        }


        @Override
        public AlloppNode getChild(int ch) {
            return ch==0 ? mlnodes[lft] : mlnodes[rgt];
        }


        @Override
        public AlloppNode getAnc() {
            return mlnodes[anc];
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
            int newch = ((MulLabNode)newchild).nodeNumber;
            if (ch == 0) {
                lft = newch;
            } else {
                rgt = newch;
            }
        }


        @Override
        public void setAnc(AlloppNode anc) {
            this.anc = ((MulLabNode)anc).nodeNumber;
        }

        @Override
        public Taxon getTaxon() {
            return taxon;
        }


        @Override
        public void setTaxon(String name) {
            taxon = new Taxon(name);
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
            lft = ((MulLabNode)c0).nodeNumber;
            mlnodes[lft].anc = nodeNumber;
            rgt = ((MulLabNode)c1).nodeNumber;
            mlnodes[rgt].anc = nodeNumber;
        }


        @Override
        public int getNumber() {
            return nodeNumber;
        }

        @Override
        public void setNumber(int nn) {
            nodeNumber = nn;
        }
    }


    private class SpSqUnion {
        public FixedBitSet spsqunion;
        public FixedBitSet spunion;
        public SpSqUnion(FixedBitSet spsqunion) {
            this.spsqunion = spsqunion;
            spunion = apsp.spsqunion2spunion(spsqunion);
        }
    }


    private class PopulationAndLineages {
        public double t[];
        public double tippop;
        public double rootpop;
        public int tipnlin;

        public PopulationAndLineages(double[] t2, double tippop, double rootpop,
                                     int tipnlin) {
            this.t = t2;
            this.tippop = tippop;
            this.rootpop = rootpop;
            this.tipnlin = tipnlin;
        }

        public double populationAt(double x) {
            final double begt = t[0];
            final double endt = t[t.length - 1];
            return ((endt-x)*tippop + (x-begt)*rootpop) / (endt-begt);
        }
    }




    /*
      * This constructor makes a single multiply labelled tree from the diploid
      * history and set of tetraploid AlloppLeggedTrees which is passed to it. It is called directly
      * by testing code.
      */
    AlloppMulLabTree(AlloppDiploidHistory adhist, ArrayList<AlloppLeggedTree> tettrees,
                     AlloppSpeciesBindings apsp, Parameter tippopvals, Parameter rootpopvals, double [] hybpopvals) {
        this.apsp = apsp;
        this.tippopvals = tippopvals;
        this.rootpopvals = rootpopvals;
        this.hybpopvals = hybpopvals;
        nofhybpopvals = tettrees.size();
        // Count tips for the tree to be made
        int ntips = 0;
        ntips += adhist.getDiploidTipCount();
        for (AlloppLeggedTree ttree : tettrees) {
            ntips += 2 * ttree.getExternalNodeCount();
        }
        // Make array of dud nodes
        mlnodes = new MulLabNode[2 * ntips - 1];
        for (int i = 0; i < mlnodes.length; i++) {
            mlnodes[i] = new MulLabNode(i);
            mlnodes[i].tetraancestor = true;
        }
        // get dip hist and tetratrees into array; unions are filled at tips
        nextn = 0;
        nextn = subtree2MulLabNodes(adhist, adhist.getRootIndex(), tettrees, apsp);
        rootn = nextn - 1;
        // fill unions from tips
        mlnodes[rootn].fillinUnionsInSubtree(apsp.numberOfSpSeqs());
        // The tetraroot flags were initialised false, set to final values by subtree2MulLabNodes()
        // The tetraancestor flags were initialised true, diploid tips later set false by subtree2MulLabNodes()
        // The intetratree flags were initialised false, tetraploid tips later set true by allopptree2MulLabNodes()
        // Now fill in tetraancestor and intetratree flags from tips.
        fillinTetraFlagsInSubtree(mlnodes[rootn]);
        makesimpletree();
    }





    // constructor for testing conversion of diploid history plus tetraploid trees to MUL-tree
    public AlloppMulLabTree(AlloppDiploidHistory adhist, ArrayList<AlloppLeggedTree> tettrees, AlloppSpeciesBindings apsp,
                            Parameter testtippopvalues, Parameter testrootpopvalues, double [] testhybpopvalues, int testcase) {
        this(adhist, tettrees, apsp, testtippopvalues, testrootpopvalues, testhybpopvalues);
        assert testcase >= 1;
        assert testcase <= 5;
        fillinpopvals();
    }


    // constructor for testing likelihood calculations.
    // Makes a particular multree with nlineages, coalheights so test can call
    // geneTreeInMULTreeLogLikelihood()
    public AlloppMulLabTree(AlloppSpeciesBindings apsp) {
        this.apsp = apsp;
        nofhybpopvals = 1;
        tippopvals = new Parameter.Default(3, .003);
        rootpopvals = new Parameter.Default(4, .001);
        hybpopvals = new double[1];
        hybpopvals[0] = .001;
        fillmlnodesforlhoodtest1();
        fillinTetraFlagsInSubtree(mlnodes[rootn]);
        mlnodes[rootn].fillinUnionsInSubtree(4);
    }


    public double testGeneTreeInMULTreeLogLikelihood() {
        return geneTreeInMULTreeLogLikelihood();
    }

    private void fillmlnodesforlhoodtest1() {
        mlnodes = new MulLabNode[7];
        for (int i = 0; i < mlnodes.length; i++) {
            mlnodes[i] = new MulLabNode(i);
        }
        // 4 tips, 2 dips, 1 tet
        mlnodes[0].taxon = new Taxon("a");
        mlnodes[1].taxon = new Taxon("b");
        mlnodes[2].taxon = new Taxon("z0");
        mlnodes[3].taxon = new Taxon("z1");
        mlnodes[2].tetraancestor = true;
        mlnodes[3].tetraancestor = true;
        mlnodes[2].intetratree = true;
        mlnodes[3].intetratree = true;
        mlnodes[2].tetraroot = true;
        mlnodes[3].tetraroot = true;
        mlnodes[0].union = new FixedBitSet(4); mlnodes[0].union.set(0);
        mlnodes[1].union = new FixedBitSet(4); mlnodes[1].union.set(1);
        mlnodes[2].union = new FixedBitSet(4); mlnodes[2].union.set(2);
        mlnodes[3].union = new FixedBitSet(4); mlnodes[3].union.set(3);
        // toplogy and times
        mlnodes[4].addChildren(mlnodes[0], mlnodes[2]);
        mlnodes[5].addChildren(mlnodes[1], mlnodes[3]);
        mlnodes[6].addChildren(mlnodes[4], mlnodes[5]);
        rootn = 6;
        mlnodes[0].height = 0.0;
        mlnodes[1].height = 0.0;
        mlnodes[2].height = 0.0;
        mlnodes[3].height = 0.0;
        mlnodes[4].height = 0.01;
        mlnodes[5].height = 0.02;
        mlnodes[6].height = 0.03;
        mlnodes[2].hybridheight = 0.005;
        mlnodes[3].hybridheight = 0.005;
        // nlineages and coalheights
        mlnodes[0].nlineages = 1;
        mlnodes[1].nlineages = 2;
        mlnodes[2].nlineages = 1;
        mlnodes[3].nlineages = 1;
        mlnodes[4].nlineages = 2;
        mlnodes[4].coalheights.add(0.015);
        mlnodes[5].nlineages = 3;
        mlnodes[5].coalheights.add(0.025);
        mlnodes[6].nlineages = 3;
        mlnodes[6].coalheights.add(0.035);
        mlnodes[6].coalheights.add(0.045);
    }



    // For testing.
    boolean mullabtreeOK() {
        int nroots = 0;
        for (int i = 0; i < mlnodes.length; i++) {
            if (mlnodes[i].anc < 0) {
                nroots++;
            }
        }
        if (nroots != 1) {
            return false;
        }
        for (int i = 0; i < mlnodes.length; i++) {
            int nparents = 0;
            for (int j = 0; j < mlnodes.length; j++) {
                if (mlnodes[j].lft == i) { nparents++; }
                if (mlnodes[j].rgt == i) { nparents++; }
            }
            if (mlnodes[i].anc < 0  &&  nparents != 0) {
                return false;
            }
            if (mlnodes[i].anc >= 0  &&  nparents != 1) {
                return false;
            }
        }
        for (int i = 0; i < mlnodes.length; i++) {
            if (mlnodes[i].getNumber() != i) {
                return false;
            }
        }
        for (int i = 0; i < mlnodes.length; i++) {
            if (mlnodes[i].lft >= 0) {
                if (mlnodes[i].rgt < 0) {
                    return false;
                }
                int lft = mlnodes[i].lft;
                int rgt = mlnodes[i].rgt;
                if (mlnodes[lft].anc != i) {
                    return false;
                }
                if (mlnodes[rgt].anc != i) {
                    return false;
                }
                if (mlnodes[i].height <= mlnodes[lft].height) {
                    return false;
                }
                if (mlnodes[i].height <= mlnodes[rgt].height) {
                    return false;
                }
            } else {
                if (mlnodes[i].height != 0) {
                    return false;
                }
            }
        }
        if (mlnodes[rootn].anc >= 0) {
            return false;
        }
        return true;
    }



    String mullabTreeAsNewick() {
        String s = TreeUtils.uniqueNewick(simptree, simptree.getRoot());
        return s;
    }

    String asText() {
        String header = "MUL-tree              height                          union                               []  tippop  []  hybpop  [] rootpop  tetroot   hybhgt nlin coalheights" + System.getProperty("line.separator");

        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + AlloppNode.Abstract.subtreeAsText(mlnodes[rootn], s, x, 0, "");
    }


    void clearCoalescences() {
        clearSubtreeCoalescences(mlnodes[rootn]);
    }

    void recordLineageCounts() {
        recordSubtreeLineageCounts(mlnodes[rootn]);
    }

    boolean coalescenceIsCompatible(double height, FixedBitSet union) {
        MulLabNode node = (MulLabNode) mlnodes[rootn].nodeOfUnionInSubtree(union);
        return (node.height <= height);
    }


    void recordCoalescence(double height, FixedBitSet union) {
        MulLabNode node = (MulLabNode) mlnodes[rootn].nodeOfUnionInSubtree(union);
        assert (node.height <= height);
        while (node.anc >= 0  &&  mlnodes[node.anc].height <= height) {
            node = mlnodes[node.anc];
        }
        node.coalheights.add(height);
    }

    void sortCoalescences() {
        for (MulLabNode node : mlnodes) {
            Collections.sort(node.coalheights);
        }
    }


    double geneTreeInMULTreeLogLikelihood() {
        fillinpopvals();
        //System.out.println(asText());
        return geneTreeInMULSubtreeLogLikelihood(mlnodes[rootn]);
    }





/*
*
* ***************************************************
* 	                Private
*/



    private void fillinTetraFlagsInSubtree(AlloppNode node) {
        if (node.nofChildren() == 2) {
            MulLabNode mnode = (MulLabNode)node;
            MulLabNode ch0 = (MulLabNode)node.getChild(0);
            MulLabNode ch1 = (MulLabNode)node.getChild(1);
            fillinTetraFlagsInSubtree(node.getChild(0));
            fillinTetraFlagsInSubtree(node.getChild(1));
            mnode.tetraancestor = (ch0.tetraancestor && ch1.tetraancestor);
            mnode.intetratree = (ch0.intetratree && !ch0.tetraroot && ch1.intetratree && !ch1.tetraroot);
        }
    }


    private int subtree2MulLabNodes(AlloppDiploidHistory adhist, int dhni, ArrayList<AlloppLeggedTree> tettrees, AlloppSpeciesBindings apsp) {
        if (adhist.getLftFromIndex(dhni) < 0) {
            int tt = adhist.getNodeTettree(dhni);
            if (tt < 0) {
                mlnodes[nextn].setTaxon(adhist.getTaxonFromIndex(dhni).getId());
                mlnodes[nextn].setUnion(apsp.taxonseqToTipUnion(mlnodes[nextn].taxon, 0));
                mlnodes[nextn].setHeight(adhist.getHeightFromIndex(dhni));
                mlnodes[nextn].tetraancestor = false;
                nextn++;
            } else {
                AlloppNode troot = (AlloppNode)tettrees.get(tt).getSlidableRoot();
                int seq = (adhist.getNodeLeg(dhni) == AlloppDiploidHistory.LegLorR.left) ? 0 : 1;
                nextn = allopptree2MulLabNodes(apsp, troot, seq);
                mlnodes[nextn-1].hybridheight = adhist.getHeightFromIndex(dhni);
                mlnodes[nextn-1].tetraroot = true;
                mlnodes[nextn-1].ttreeindex = tt;
            }
        } else {
            nextn = subtree2MulLabNodes(adhist, adhist.getLftFromIndex(dhni), tettrees, apsp);
            int c0 = nextn - 1;
            nextn = subtree2MulLabNodes(adhist, adhist.getRgtFromIndex(dhni), tettrees, apsp);
            int c1 = nextn - 1;
            mlnodes[nextn].addChildren(mlnodes[c0], mlnodes[c1]);
            mlnodes[nextn].setHeight(adhist.getHeightFromIndex(dhni));
            nextn++;
        }
        return nextn;
    }



    private int allopptree2MulLabNodes(AlloppSpeciesBindings apsp,
                                       AlloppNode snode, int seq) {
        if (snode.nofChildren() == 0) {
            mlnodes[nextn].setTaxon(snode.getTaxon().getId() + seq);
            mlnodes[nextn].setUnion(apsp.taxonseqToTipUnion(snode.getTaxon(), seq));
            mlnodes[nextn].intetratree = true;
        } else {
            nextn = allopptree2MulLabNodes(apsp, snode.getChild(0), seq);
            int c0 = nextn - 1;
            nextn = allopptree2MulLabNodes(apsp, snode.getChild(1), seq);
            int c1 = nextn - 1;
            mlnodes[nextn].addChildren(mlnodes[c0], mlnodes[c1]);
        }
        mlnodes[nextn].setHeight(snode.getHeight());
        nextn++;
        return nextn;
    }





    private void makesimpletree() {
        SimpleNode[] snodes = new SimpleNode[mlnodes.length];
        for (int n = 0; n < mlnodes.length; n++) {
            snodes[n] = new SimpleNode();
        }
        makesimplesubtree(snodes, 0, mlnodes[rootn]);
        simptree = new SimpleTree(snodes[mlnodes.length-1]);
    }


    private int makesimplesubtree(SimpleNode[] snodes, int nextsn, MulLabNode mnode) {
        if (mnode.lft < 0) {
            snodes[nextsn].setTaxon(new Taxon(mnode.taxon.getId()));
        } else {
            nextsn = makesimplesubtree(snodes, nextsn, mlnodes[mnode.lft]);
            int subtree0 = nextsn-1;
            nextsn = makesimplesubtree(snodes, nextsn, mlnodes[mnode.rgt]);
            int subtree1 = nextsn-1;
            snodes[nextsn].addChild(snodes[subtree0]);
            snodes[nextsn].addChild(snodes[subtree1]);
        }
        snodes[nextsn].setHeight(mnode.height);
        String tti;
        if (mnode.ttreeindex < 0) {
            tti = "X";
        } else {
            tti = "T" + mnode.ttreeindex;
        }
        snodes[nextsn].setAttribute("tti", tti);
        if (mnode.hybridheight >= 0.0) {
            snodes[nextsn].setAttribute("hybhgt", mnode.hybridheight);
        }
        return nextsn+1;
    }


    private void clearSubtreeCoalescences(MulLabNode node) {
        if (node.lft >= 0) {
            clearSubtreeCoalescences(mlnodes[node.lft]);
            clearSubtreeCoalescences(mlnodes[node.rgt]);
        }
        node.coalheights.clear();
    }


    private void recordSubtreeLineageCounts(MulLabNode node) {
        if (node.lft < 0) {
            node.nlineages = apsp.nLineages(apsp.spseqindex2sp(union2spseqindex(node.union)));
        } else {
            node.nlineages = 0;
            recordSubtreeLineageCounts(mlnodes[node.lft]);
            recordSubtreeLineageCounts(mlnodes[node.rgt]);
            node.nlineages += mlnodes[node.lft].nlineages - mlnodes[node.lft].coalheights.size();
            node.nlineages += mlnodes[node.rgt].nlineages - mlnodes[node.rgt].coalheights.size();
        }
    }



    /*
      * This copies population values in the Parameter popvalues
      * to nodes in the AlloppMulLabTree. The population values are
      * per-species-clade (per-branch in network), but of course more than
      * one node in AlloppMulLabTree may correspond to the same species.
      *
      * The other complications are that tips are different from internal
      * nodes, and that nodes which roots of tetratrees or just below,
      * as well as the root are special cases.
      *
      * It collects unions (which represent sets whose elements
      * identify a species and a sequence) from the nodes and then
      * sorts them primarily using identities of the species, so
      * that sets of node with same species clade are grouped together. The sort
      * also puts the node sets corresponding to tips first in the array and sorts
      * nodes within node sets in a well-defined way.
      * This mainly does what is required, since nodes with the same
      * species clade are treated the same.
      *
      * fillinpopvalsforspunion() deals with a set of nodes
      * with same species clade.
      */
    private void fillinpopvals() {
        ArrayList<SpSqUnion> unionarraylist = new ArrayList<SpSqUnion>();
        for (int n = 0; n < mlnodes.length; n++) {
            unionarraylist.add(new SpSqUnion(mlnodes[n].union));
        }
        Collections.sort(unionarraylist, SPUNION_ORDER);
        SpSqUnion[] unionarray = new SpSqUnion[unionarraylist.size()];
        unionarray = unionarraylist.toArray(unionarray);
        PopValIndices pvis = new PopValIndices(0,0,0);
        // set all pop indices to dud values
        for (int n = 0; n < mlnodes.length; n++) {
            mlnodes[n].tippopindex = -1;
            mlnodes[n].rootpopindex = -1;
            mlnodes[n].hybpopindex = -1;
        }
        int noftippopvals = tippopvals.getDimension();
        int nofrootpopvals = rootpopvals.getDimension();
        for (int n0 = 0; n0 < unionarray.length; ) {
            int n1 = n0+1;
            while (n1 < unionarray.length && unionarray[n1].spunion.equals(unionarray[n0].spunion)) {
                n1++;
            }
            assert pvis.tipp <= noftippopvals;
            if (!(pvis.rootp <= nofrootpopvals)) {
                System.out.println("BUG in fillinpopvals()");
            }
            assert pvis.rootp <= nofrootpopvals;
            assert pvis.hybp <= nofhybpopvals;
            pvis = fillinpopvalsforspunion(unionarray, n0, n1, pvis);
            n0 = n1;
        }
        if (pvis.tipp != noftippopvals  || pvis.rootp != nofrootpopvals  ||  pvis.hybp != nofhybpopvals) {
            System.out.println("BUG in fillinpopvals()");
        }
        assert pvis.tipp == noftippopvals;
        assert pvis.rootp == nofrootpopvals;
        assert pvis.hybp == nofhybpopvals;
    }


    private class PopValIndices {
        public int tipp;
        public int rootp;
        public int hybp;
        PopValIndices(int tipp, int rootp, int hybp) {
            this.tipp = tipp;
            this.rootp = rootp;
            this.hybp = hybp;
        }
    }


    private PopValIndices fillinpopvalsforspunion(SpSqUnion[] unionarray, int n0, int n1, PopValIndices pvis) {
        int n = n1-n0;
        MulLabNode nodeset[] = new MulLabNode[n];
        // get set of nodes with same species clade
        for (int i = n0; i < n1; i++) {
            nodeset[i-n0] = (MulLabNode) mlnodes[rootn].nodeOfUnionInSubtree(unionarray[i].spsqunion);
        }
        // Hybridization pops. Nodes in pairs, shared pop.
        if (nodeset[0].tetraroot) {
            assert n >= 2;
            assert nodeset[1].tetraroot;
            for (int i = 0; i < 2; i++) {
                nodeset[i].hybpopindex = pvis.hybp;
            }
            pvis.hybp++;
        }
        // Tip pops. Either dip or tet; in latter case shared pop.
        int ntips = 0;
        for (int i = 0; i < n; i++) {
            if (nodeset[i].lft < 0) {
                nodeset[i].tippopindex = pvis.tipp;
                ntips++;
            }
        }
        assert ntips <= 2;
        if (ntips > 0) {
            pvis.tipp++;
        }

        // Root pops within tetra trees. Nodes in pairs.
        if (nodeset[0].intetratree  &&  !nodeset[0].tetraroot) {
            if (n != 2) {
                System.out.println("BUG in fillinpopvalsforspunion() 2");
            }
            assert n == 2;
            assert (nodeset[1].intetratree  &&  !nodeset[1].tetraroot);
            for (int i = 0; i < 2; i++) {
                nodeset[i].rootpopindex = pvis.rootp;
            }
            pvis.rootp++;
        } else {
            // Root pops NOT within tetra trees. Here the nodeset[] structure is not useful.
            // Here, pop is shared with sibling if either is ancestor to tetraploids only
            for (int i = 0; i < n; i++) {
                MulLabNode node = nodeset[i];
                if (node.anc >= 0) {
                    MulLabNode sibling = siblingOfNode(node);
                    if (node.tetraancestor  ||  sibling.tetraancestor) {
                        if (sibling.rootpopindex >= 0) {
                            node.rootpopindex = sibling.rootpopindex;
                        } else {
                            node.rootpopindex = pvis.rootp;
                            pvis.rootp++;
                        }
                    } else {
                        node.rootpopindex = pvis.rootp;
                        pvis.rootp++;
                    }
                }
            }
        }

        return pvis;
    }



    private MulLabNode siblingOfNode(MulLabNode node) {
        assert node.anc >= 0;
        MulLabNode sibling;
        if (mlnodes[ mlnodes[node.anc].lft ] == node) {
            sibling = mlnodes[ mlnodes[node.anc].rgt ];
        } else {
            assert mlnodes[ mlnodes[node.anc].rgt ] == node;
            sibling = mlnodes[ mlnodes[node.anc].lft ];
        }
        return sibling;
    }



    /*
      * Visits each node in MULtree and accumulates LogLikelihood
      */
    private double geneTreeInMULSubtreeLogLikelihood(MulLabNode node) {
        double loglike = 0.0;
        if (node.lft >= 0) {
            loglike += geneTreeInMULSubtreeLogLikelihood(mlnodes[node.lft]);
            loglike += geneTreeInMULSubtreeLogLikelihood(mlnodes[node.rgt]);
        }
        loglike += branchLLInMULtree(node);

        return loglike;
    }


    /*
      * Does likelihood calculation for a single node in the case
      * of two diploids.
      */
    private double branchLLInMULtree(MulLabNode node) {
        double loglike = 0.0;

        double tippop = 0.0;
        if (node.lft < 0) {
            tippop = node.tippop();
        } else {
            tippop = mlnodes[node.lft].rootpop() + mlnodes[node.rgt].rootpop();
        }

        PopulationAndLineages pal;
        double t[];
        if (node.tetraroot) {
            // since hybridization
            int nsince = 0;
            for ( ; nsince < node.coalheights.size() &&
                    node.coalheights.get(nsince) < node.hybridheight; nsince++) {}
            t = new double[nsince + 2];
            t[0] = node.height;
            for (int i = 0; i < nsince; i++) {
                t[i+1] = node.coalheights.get(i);
            }
            t[t.length-1] = node.hybridheight;
            pal = new PopulationAndLineages(t, tippop, node.hybpop(), node.nlineages);
            loglike += limbLogLike(pal);
            if (Double.isNaN(loglike)) {
                System.out.println("BUG in branchLLInMULtree");
            }
            // before hybridization
            int nbefore = node.coalheights.size() - nsince;
            t = new double[nbefore + 2];
            t[0] = node.hybridheight;
            for (int i = 0; i < nbefore; i++) {
                t[i+1] = node.coalheights.get(nsince+i);
            }
            t[t.length-1] = mlnodes[node.anc].height;

            pal = new PopulationAndLineages(t, node.rootpop(), node.rootpop(), node.nlineages - nsince);
            loglike += limbLogLike(pal);
        } else if (node.anc < 0) {
            t = new double[node.coalheights.size() + 2];
            t[0] = node.height;
            t[t.length-1] = apsp.maxGeneTreeHeight();
            for (int i = 0; i < node.coalheights.size(); i++) {
                t[i+1] = node.coalheights.get(i);
            }
            pal = new PopulationAndLineages(t, tippop, tippop, node.nlineages);
            loglike += limbLogLike(pal);
        } else {
            t = new double[node.coalheights.size() + 2];
            t[0] = node.height;
            t[t.length-1] = mlnodes[node.anc].height;
            for (int i = 0; i < node.coalheights.size(); i++) {
                t[i+1] = node.coalheights.get(i);
            }
            pal = new PopulationAndLineages(t, tippop, node.rootpop(), node.nlineages);
            loglike += limbLogLike(pal);
        }
        if (Double.isNaN(loglike)) {
            System.out.println("BUG in branchLLInMULtree");
        }
        return loglike;
    }




    /*
      * limbLogLike calculates the log-likelihood for
      * the coalescences at t[1],t[2],...t[k] within a limb
      * from t[0] to t[k+1]. ('limb' means a branch or part of one.)
      */
    private double limbLogLike(PopulationAndLineages pal) {
        double loglike = 0.0;
        int k = pal.t.length - 2;
        for (int i = 1; i <= k; i++) {
            loglike -= Math.log(pal.populationAt(pal.t[i]));
        }
        for (int i = 0; i <= k; i++) {
            final double y = (pal.tipnlin-i) * (pal.tipnlin-i-1) / 2 ;
            final double z = limbLinPopIntegral(pal, pal.t[i], pal.t[i+1]);
            loglike -= y * z;
        }
        if (Double.isNaN(loglike)) {
            System.out.println("BUG in limbLogLike");
        }
        return loglike;
    }



    // integral from t0 to t1 of (endt-begt)/((endt-x)begPop + (x-begt)endPop)
    // with respect to x
    private double limbLinPopIntegral(PopulationAndLineages b, double t0, double t1) {
        final double begt = b.t[0];
        final double endt = b.t[b.t.length-1];
        if (b.rootpop < 1E-20) {
            System.out.println("Underflow in limbLinPopIntegral()");
        }
        final double d = b.rootpop - b.tippop;
        final double c = endt * b.tippop - begt * b.rootpop;
        final double x = Math.abs(d / b.tippop);
        if (x > 0.001) {
            return ((endt - begt) / d) * Math.log((c + d * t1) / (c + d * t0));
        } else {
            double y = d * (t1 - t0) / (c + d * t0);
            double ys = (1.0 - y/2.0 + y*y * (1.0/3.0 - y/4.0 + y*y * (1.0/5.0 - y/6.0)));
            return ((endt - begt) * (t1 - t0) / (c + d * t0)) * ys;
        }
    }




    private static int union2spseqindex(FixedBitSet union) {
        assert union.cardinality() == 1;
        return union.nextOnBit(0);
    }




    /*
      * This is for ordering the unions in the nodes of the AlloppMulLabTree.
      * Those unions are of (species, sequence) pairs.
      *
      * The comparator sorts the unions of (species, sequence) pairs (SpSqUnions)
      * so that all unions containing the same set of species (ignoring sequence)
      * are grouped together. Call the sets of SpSqUnions for the same species a
      * `group'. There can be 1,2 or 3 SpSqUnions in a group.
      *
      * The groups are sorted in order of increasing number of species (clade size).
      * All groups for a single species (a tip in the network) come first, then
      * those groups for two species, and so on to the root for all species.
      * For groups that have equal numbers of species, a lexicographical
      * ordering using species indices is used.
      *
      * Within each group, species and sequence information is used to sort the 1 to 3
      * SpSqUnions. The size of the `clade' of (species, sequence) pairs is used
      * first in the comparison, which ensures that the three nodes with the same species
      * - corresponding to two roots of tetratrees in the AlloppMulLabTree plus a leg-join -
      * are ordered so that the two roots come first.
      *
      */

    static final Comparator<SpSqUnion> SPUNION_ORDER = new Comparator<SpSqUnion>() {
        public int compare(SpSqUnion a, SpSqUnion b) {
            int ac = a.spunion.cardinality();
            int bc = b.spunion.cardinality();
            if (ac != bc) {
                return ac-bc;
            } else {
                int an = a.spunion.nextOnBit(0);
                int bn = b.spunion.nextOnBit(0);
                while (an >= 0 || bn >= 0) {
                    if (an != bn) {
                        return an-bn;
                    }
                    an = a.spunion.nextOnBit(an+1);
                    bn = b.spunion.nextOnBit(bn+1);
                }
                // spunions compare equal; do spsqunions
                ac = a.spsqunion.cardinality();
                bc = b.spsqunion.cardinality();
                if (ac != bc) {
                    return ac-bc;
                } else {
                    an = a.spsqunion.nextOnBit(0);
                    bn = b.spsqunion.nextOnBit(0);
                    while (an >= 0 || bn >= 0) {
                        if (an != bn) {
                            return an-bn;
                        }
                        an = a.spsqunion.nextOnBit(an+1);
                        bn = b.spsqunion.nextOnBit(bn+1);
                    }
                    return 0;
                }
            }
        }
    };




}
