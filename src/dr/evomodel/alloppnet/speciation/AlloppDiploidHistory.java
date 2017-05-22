/*
 * AlloppDiploidHistory.java
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
import java.util.Formatter;
import java.util.Locale;
import java.util.Stack;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.alloppnet.tree.SlidableTree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;
import dr.evomodel.alloppnet.util.AlloppMisc;
import jebl.util.FixedBitSet;


/**
 * AlloppDiploidHistory represents part of the network before hybridizations.
 * It is basically a tree with some tips representing diploid species at present time
 * and others representing the points at which hybridization occurs.
 *
 * @author Graham Jones
 *         Date: 13/03/2012
 */


/*
 * This a tree with some tips representing diploid species at present time
 * and others representing the points at which hybridization (to form a
 * tetraploid) occurs. The latter are in pairs and have times before present,
 * and I call them `hydridization tips' or `hyb-tips'. From the point of
 * view of a tetraploid tree these hyb-tips are feet at the ends of legs.
 *
 * The purpose of this class is to represent the part of the network before
 * hybridizations (=`diploid history') in a form (the array of DipHistNode's)
 * that can be subjected to Mau-type moves.
 *
 * The diploid history is constructed from a diploid tree and a single
 * tetraploid tree initially. MCMC moves can change the number of tetraploid
 * trees.
 *
 */



public class AlloppDiploidHistory implements SlidableTree {
    private DipHistNode[] dhnodes;
    private int rootn;
    private int nextn;
    private AlloppSpeciesBindings apsp;


    /******************************** Inner classes ****************************************/

    public enum LegLorR {left, right, dud};


    // Small class for returning two heights. For ChangeNumHybs move
    public class FootAncHeights {
        public double anchgt;
        public double ancanchgt;
        FootAncHeights(double anchgt, double ancanchgt) {
            this.anchgt = anchgt;
            this.ancanchgt = ancanchgt;
        }
    }

    /*
    * parent, child[] implement the tree topology.
    *
    * height is the node height; can be > 0 for hyb-tips.
    *
    * union is a spseq-union. At diploid tips they are used normally.
    * Hyb-tips have unions derived from the tetratree and leg index (0 or 1).
    * Unions are taken towards the root.
    *
    * For a hyb-tip, tettree specifies the index of the tetratree whose
    * root comes from it. It is not used for other nodes.
    */
    private class DipHistNode extends AlloppNode.Abstract implements AlloppNode, NodeRef {
        private int anc;
        private int lft;
        private int rgt;
        private double height;
        private Taxon taxon;
        private FixedBitSet union;
        private int tettree;
        private LegLorR leg;
        private int nodeNumber;

        // dud constuctor
        DipHistNode(int nn) {
            anc = -1;
            lft = -1;
            rgt = -1;
            height = -1.0;
            taxon = new Taxon("");
            union = null;
            tettree = -1;
            leg = LegLorR.dud;
            nodeNumber = nn;
        }

        // copy constructor
        public DipHistNode(DipHistNode node) {
            anc = node.anc;
            lft = node.lft;
            rgt = node.rgt;
            nodeNumber = node.nodeNumber;
            copyNonTopologyFields(node);
        }


        // no-topology constructor. Copies all fields of argument node,
        // except the anc, lft, rgt fields which are set to `unknown',
        // and the index field, which is set to nn
        public DipHistNode(int nn, DipHistNode node) {
            anc = -1;
            lft = -1;
            rgt = -1;
            nodeNumber = nn;
            copyNonTopologyFields(node);
        }

        private void copyNonTopologyFields(DipHistNode node) {
            height = node.height;
            taxon = new Taxon(node.taxon.getId());
            if (node.union == null) {
                union = null;
            } else {
                union = new FixedBitSet(node.union);
            }
            tettree = node.tettree;
            leg = node.leg;
        }


        @Override
        public AlloppNode getChild(int ch) {
            return ch==0 ? dhnodes[lft] : dhnodes[rgt];
        }

        @Override
        public AlloppNode getAnc() {
            return dhnodes[anc];
        }


        @Override
        public double getHeight() {
            return height;
        }

        @Override
        public Taxon getTaxon() {
            return taxon;
        }

        @Override
        public FixedBitSet getUnion() {
            return union;
        }

        @Override
        public void setChild(int ch, AlloppNode newchild) {
            int newch = ((DipHistNode)newchild).nodeNumber;
            if (ch == 0) {
                lft = newch;
            } else {
                rgt = newch;
            }
        }

        @Override
        public void setAnc(AlloppNode anc) {
            this.anc = ((DipHistNode)anc).nodeNumber;
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
            lft = ((DipHistNode)c0).nodeNumber;
            dhnodes[lft].anc = nodeNumber;
            rgt = ((DipHistNode)c1).nodeNumber;
            dhnodes[rgt].anc = nodeNumber;
        }

        @Override
        public String asText(int indentlen) {
            StringBuilder s = new StringBuilder();
            Formatter formatter = new Formatter(s, Locale.US);
            if (lft < 0) {
                String nodename;
                if (tettree >= 0) {
                    nodename = String.valueOf(tettree);
                    String legtext = "*";
                    if (leg == LegLorR.left) {
                        legtext = "L";
                    }
                    if (leg == LegLorR.right) {
                        legtext = "R";
                    }
                    nodename += legtext;
                } else {
                    nodename = taxon.getId();
                }
                formatter.format("%s ", nodename);
            } else {
                formatter.format("%s ", "+");
            }
            while (s.length() < 25-indentlen) {
                formatter.format("%s", " ");
            }
            formatter.format("%s ", AlloppMisc.nonnegIn8Chars(height));
            return s.toString();
        }


        @Override
        public int nofChildren() {
            return (lft < 0) ? 0 : 2;
        }

        @Override
        public int getNumber() {
            return nodeNumber;
        }

        @Override
        public void setNumber(int n) {
            nodeNumber = n;

        }



    }


    /******************************** Constructors ****************************************/



    private class JoiningNode {
        int nn;
        boolean hasdip;
        JoiningNode(int nn, boolean hasdip) {
            this.nn = nn;
            this.hasdip = hasdip;
        }
    }


    /*
      * Constructor makes initial random AlloppDiploidHistory from a diploid tree
      * and previously constructed tetratrees.
      * apsp needed for speciesseqEmptyUnion(). grjtodo-oneday really needed?
      *
      */
    AlloppDiploidHistory(Taxon [] dipspp, ArrayList<AlloppLeggedTree> tettrees,
                         boolean diploidrootisroot, double rate, AlloppSpeciesBindings apsp) {
        this.apsp = apsp;
        int nofdips = dipspp.length;
        int ntips = nofdips + 2 * tettrees.size();
        // Make array of dud nodes
        dhnodes = new DipHistNode[2 * ntips - 1];
        for (int i = 0; i < dhnodes.length; i++) {
            dhnodes[i] = new DipHistNode(i);
        }
        // Make tips for diploids and add them to tojoin[]
        ArrayList<JoiningNode> tojoin = new ArrayList<JoiningNode>();
        for (nextn = 0; nextn < nofdips; nextn++) {
            dhnodes[nextn].taxon = new Taxon(dipspp[nextn].getId());
            dhnodes[nextn].height = 0.0;
            tojoin.add(new JoiningNode(nextn, true));
        }
        // Make two tips for each tet tree, and add to tojoin[]
        for (int tt = 0; tt < tettrees.size(); tt++) {
            AlloppLeggedTree ttree = tettrees.get(tt);
            double tettipsheight = ttree.getRootHeight() + MathUtils.nextExponential(rate);
            dhnodes[nextn].tettree = tt;
            dhnodes[nextn].leg = LegLorR.left;
            ttree.setDiphistLftLeg(nextn);
            dhnodes[nextn].height = tettipsheight;
            tojoin.add(new JoiningNode(nextn, false));
            nextn++;
            dhnodes[nextn].tettree = tt;
            dhnodes[nextn].leg = LegLorR.right;
            ttree.setDiphistRgtLeg(nextn);
            dhnodes[nextn].height = tettipsheight;
            tojoin.add(new JoiningNode(nextn, false));
            nextn ++;
        }

        // build the tree from tips, taking care not to join the last diploids before all tets used up
        for (int i = 0; i < ntips-1; i++) {
            int numtojoin = tojoin.size();
            int numwithdips = 0;
            for (JoiningNode jn : tojoin) {
                numwithdips += jn.hasdip ? 1 : 0;
            }
            int j = MathUtils.nextInt(numtojoin);
            JoiningNode child0 = tojoin.get(j);
            // if diploidrootisroot, and down to 2 diploids, and still a tet, ensure ch0 has no dips
            if (diploidrootisroot  &&  numtojoin > 2  &&  numwithdips == 2) {
                while (child0.hasdip) {
                    j = MathUtils.nextInt(numtojoin);
                    child0 = tojoin.get(j);
                }
            }
            tojoin.remove(j);
            int k = MathUtils.nextInt(numtojoin-1);
            JoiningNode child1 = tojoin.get(k);
            tojoin.remove(k);
            dhnodes[nextn].lft = child0.nn;   dhnodes[child0.nn].anc = nextn;
            dhnodes[nextn].rgt = child1.nn;   dhnodes[child1.nn].anc = nextn;
            double maxchhgt = Math.max(dhnodes[child0.nn].height, dhnodes[child1.nn].height);
            dhnodes[nextn].height = maxchhgt + MathUtils.nextExponential(numtojoin*rate);
            tojoin.add(new JoiningNode(nextn, child0.hasdip || child1.hasdip));
            nextn++;
        }
        rootn = nextn - 1;
        assert diphistOK(diploidrootisroot);
        makesimpletree();
    }




    // copy constructor.
    public AlloppDiploidHistory(AlloppDiploidHistory x) {
        dhnodes = new DipHistNode[x.dhnodes.length];
        for (int n = 0; n < dhnodes.length; n++) {
            dhnodes[n] = new DipHistNode(x.dhnodes[n]);
        }
        rootn = x.rootn;
        nextn = x.nextn;
        apsp = x.apsp;
    }


    //  constructor for testing.
    public AlloppDiploidHistory(SimpleNode[] snodes, int sroot, ArrayList<AlloppLeggedTree> tettrees,
                                boolean diprootisroot, AlloppSpeciesBindings apsp) {
        this.apsp = apsp;
        // Make array of dud nodes
        dhnodes = new DipHistNode[snodes.length];
        for (int i = 0; i < dhnodes.length; i++) {
            dhnodes[i] = new DipHistNode(i);
        }
        // Convert snodes tree into dhnodes tree
        nextn = 0;
        simpletree2dhtesttree(snodes[sroot]);
        rootn = nextn - 1;

        // link the hyb-tips to tet roots.
        for (int n = 0; n < dhnodes.length; n++) {
            String tipname = dhnodes[n].taxon.getId();
            LegLorR leg = LegLorR.dud;
            int tt = -1;
            if (tipname.contains("L")) { leg = LegLorR.left; }
            if (tipname.contains("R")) { leg = LegLorR.right; }
            if (tipname.contains("0")) { tt = 0; }
            if (tipname.contains("1")) { tt = 1; }
            if (tipname.contains("2")) { tt = 2; }
            if (leg != LegLorR.dud) {
                assert tt >= 0;
                dhnodes[n].tettree = tt;
                dhnodes[n].leg = leg;
                if (leg == LegLorR.left) { tettrees.get(tt).setDiphistLftLeg(n); }
                if (leg == LegLorR.right) { tettrees.get(tt).setDiphistRgtLeg(n); }
            }
        }
        dhnodes[rootn].fillinUnionsInSubtree(apsp.numberOfSpSeqs());
        assert diphistOK(diprootisroot);
        makesimpletree();
    }



    /******************************** SlidableTree implementation ****************************************/


    @Override
    public NodeRef getSlidableRoot() {
        assert dhnodes[rootn].anc < 0;
        return dhnodes[rootn];
    }


    @Override
    public void replaceSlidableRoot(NodeRef root) {
        rootn = root.getNumber();
        dhnodes[rootn].anc = -1;
    }


    @Override
    public int getSlidableNodeCount() {
        return dhnodes.length;
    }


    @Override
    public double getSlidableNodeHeight(NodeRef node) {
        return dhnodes[node.getNumber()].getHeight();
    }

    @Override
    public Taxon getSlidableNodeTaxon(NodeRef node) {
        return dhnodes[node.getNumber()].getTaxon();
    }

    @Override
    public void setSlidableNodeHeight(NodeRef node, double height) {
        dhnodes[node.getNumber()].height = height;
    }


    @Override
    public boolean isExternalSlidable(NodeRef node) {
        return (dhnodes[node.getNumber()].lft < 0);
    }



    @Override
    public NodeRef getSlidableChild(NodeRef node, int j) {
        int n = node.getNumber();
        return j == 0 ? dhnodes[ dhnodes[n].lft ] : dhnodes[ dhnodes[n].rgt ];
    }



    @Override
    public void replaceSlidableChildren(NodeRef node, NodeRef lft, NodeRef rgt) {
        int nn = node.getNumber();
        int lftn = lft.getNumber();
        int rgtn = rgt.getNumber();
        assert dhnodes[nn].lft >= 0;
        dhnodes[nn].lft = lftn;
        dhnodes[nn].rgt = rgtn;
        dhnodes[lftn].anc = dhnodes[nn].nodeNumber;
        dhnodes[rgtn].anc = dhnodes[nn].nodeNumber;
    }



    String asText() {
        String header = "Diploid history            height" + System.getProperty("line.separator");
        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + AlloppNode.Abstract.subtreeAsText(dhnodes[rootn], s, x, 0, "");
    }


    int getInternalNodeCount() {
        return (dhnodes.length - 1) / 2;
    }


    int getDiploidTipCount() {
        int ndiptips = 0;
        for (int i = 0; i < dhnodes.length; i++) {
            if (dhnodes[i].lft < 0  &&  dhnodes[i].tettree < 0) {
                ndiptips++;
            }
        }
        return ndiptips;
    }


    public ArrayList<Integer> collectFeet() {
        ArrayList<Integer> feet = new ArrayList<Integer>();
        for (int i = 0; i < dhnodes.length; i++) {
            if (dhnodes[i].tettree >= 0) {
                feet.add(i);
            }
        }
        return feet;
    }


    public boolean tipIsDiploidTip(NodeRef node) {
        assert dhnodes[node.getNumber()].lft < 0;
        return (dhnodes[node.getNumber()].tettree < 0);
    }

    // For move that merges two tettrees. This is part of
    // test that they can be merged.
    public boolean tettreesShareLegs(AlloppLeggedTree ttree1, AlloppLeggedTree ttree2) {
        int lftleg1 = ttree1.getDiphistLftLeg();
        int lftleg1anc = dhnodes[lftleg1].anc;
        int lftleg2 = ttree2.getDiphistLftLeg();
        int lftleg2anc = dhnodes[lftleg2].anc;
        int rgtleg1 = ttree1.getDiphistRgtLeg();
        int rgtleg1anc = dhnodes[rgtleg1].anc;
        int rgtleg2 = ttree2.getDiphistRgtLeg();
        int rgtleg2anc = dhnodes[rgtleg2].anc;
        boolean llrr = ((lftleg1anc == lftleg2anc)  &&  (rgtleg1anc == rgtleg2anc));
        return llrr;
        /*
        boolean lrrl = ((lftleg1anc == rgtleg2anc)  &&  (rgtleg1anc == lftleg2anc));
        return  (llrr || lrrl);  */
    }


    // For (old?) move that merges two tettrees. This is for Hastings ratio calculation.
    public double intervalOfFoot(AlloppLeggedTree ttree, boolean left) {
        double hybh = getHybHeight(ttree);
        int foot = left ? ttree.getDiphistLftLeg() : ttree.getDiphistRgtLeg();
        int footanc = dhnodes[foot].anc;
        assert footanc >= 0;
        int footancanc = dhnodes[footanc].anc;
        assert footancanc >= 0;
        return (dhnodes[footancanc].height  -  hybh);
    }


    // For move that merges two tettrees. This is for Hastings ratio calculation.
    public FootAncHeights intervalOfFootAncestor(AlloppLeggedTree ttree, LegLorR leg) {
        double hybh = getHybHeight(ttree);
        int foot = (leg==LegLorR.left) ? ttree.getDiphistLftLeg() : ttree.getDiphistRgtLeg();
        assert hybh == dhnodes[foot].height;
        int footanc = dhnodes[foot].anc;
        int footancanc = dhnodes[footanc].anc;
        assert footancanc >= 0;
        return new FootAncHeights(dhnodes[footanc].height, dhnodes[footancanc].height);
    }



    public void setHybridHeight(AlloppLeggedTree ttree, double newh) {
        int foot1 = ttree.getDiphistLftLeg();
        int foot2 = ttree.getDiphistRgtLeg();
        dhnodes[foot1].height = dhnodes[foot2].height = newh;
    }



    // For move that merges two tettrees. This is for after merge.
    // It removes two tips that were joined to ttree
    public void removeFeet(AlloppSpeciesNetworkModel apspnet, AlloppLeggedTree ttree) {
        DipHistNode [] tmpnodes  = new DipHistNode[dhnodes.length];
        for (int i = 0;  i < tmpnodes.length; i++) {
            tmpnodes[i] = new DipHistNode(dhnodes[i]);
        }
        removeTip(ttree.getDiphistLftLeg(), tmpnodes);
        removeTip(ttree.getDiphistRgtLeg(), tmpnodes);
        dhnodes = new DipHistNode[tmpnodes.length - 4];
        nextn = 0;
        buildSubtreeFromNodes(apspnet, tmpnodes, rootn);
        rootn = nextn - 1;
    }


    // Adds two hyb-tips for tettree tt1. these are joined to legs of tettree tt2
    public void addTwoDipTips(AlloppSpeciesNetworkModel apspnet, int tt1, int tt2, double lfthgt, double rgthgt, double tiphgt) {
        AlloppLeggedTree tettree1 = apspnet.getTetraploidTree(tt1);
        AlloppLeggedTree tettree2 = apspnet.getTetraploidTree(tt2);
        int lftleg = tettree2.getDiphistLftLeg();
        int rgtleg = tettree2.getDiphistRgtLeg();
        int oldn = dhnodes.length;
        DipHistNode [] tmpnodes = new DipHistNode[oldn + 4];
        for (int n = 0; n < oldn; n++) {
            tmpnodes[n] = new DipHistNode(dhnodes[n]);
        }

        // two new diploid tips for tettree tt1
        tmpnodes[oldn] = new DipHistNode(oldn);
        tmpnodes[oldn].height = tiphgt;
        tmpnodes[oldn].tettree = tt1;
        tmpnodes[oldn].leg = LegLorR.left;
        tettree1.setDiphistLftLeg(oldn);

        tmpnodes[oldn+1] = new DipHistNode(oldn+1);
        tmpnodes[oldn+1].height = tiphgt;
        tmpnodes[oldn+1].tettree = tt1;
        tmpnodes[oldn+1].leg = LegLorR.right;
        tettree1.setDiphistRgtLeg(oldn+1);

        // two new nodes to go into existing branches
        tmpnodes[oldn+2] = new DipHistNode(oldn+2);
        tmpnodes[oldn+2].height = lfthgt;
        tmpnodes[oldn+2].anc = tmpnodes[lftleg].anc;
        tmpnodes[oldn+2].lft = oldn;
        tmpnodes[oldn+2].rgt = lftleg;

        tmpnodes[oldn+3] = new DipHistNode(oldn+3);
        tmpnodes[oldn+3].height = rgthgt;
        tmpnodes[oldn+3].anc = tmpnodes[rgtleg].anc;
        tmpnodes[oldn+3].lft = oldn+1;
        tmpnodes[oldn+3].rgt = rgtleg;

        // divide branch by pointing to new nodes oldn+2,+3
        int lftfootanc = tmpnodes[lftleg].anc;
        if (tmpnodes[lftfootanc].lft == lftleg) {
            tmpnodes[lftfootanc].lft = oldn+2;
        } else {
            assert tmpnodes[lftfootanc].rgt == lftleg;
            tmpnodes[lftfootanc].rgt = oldn+2;
        }

        int rgtfootanc = tmpnodes[rgtleg].anc;
        if (tmpnodes[rgtfootanc].lft == rgtleg) {
            tmpnodes[rgtfootanc].lft = oldn+3;
        } else {
            assert tmpnodes[rgtfootanc].rgt == rgtleg;
            tmpnodes[rgtfootanc].rgt = oldn+3;
        }

        dhnodes = new DipHistNode[tmpnodes.length];
        nextn = 0;
        buildSubtreeFromNodes(apspnet, tmpnodes, rootn);
        rootn = nextn - 1;
    }


    // For move that splits a tetree
    public double getAncHeight(int n) {
        int nanc = dhnodes[n].anc;
        assert nanc >= 0;
        return dhnodes[nanc].height;
    }


    // next bunch for AlloppMulLabTree constructor
    int getRootIndex() {
        return rootn;
    }

    double getHeightFromIndex(int n) {
        return dhnodes[n].height;
    }

    int getLftFromIndex(int n) {
        return dhnodes[n].lft;
    }
    int getRgtFromIndex(int n) {
        return dhnodes[n].rgt;
    }
    Taxon getTaxonFromIndex(int n) {
        return dhnodes[n].taxon;
    }


    double getRootHeight() {
        return dhnodes[rootn].height;
    }


    public double getHybHeight(AlloppLeggedTree tettree) {
        return dhnodes[tettree.getDiphistLftLeg()].height;
    }


    // for prior lhood
    void collectInternalAndHybHeights(ArrayList<Double> heights) {
        for (DipHistNode node : dhnodes) {
            if (node.lft >= 0) {
                heights.add(node.height);
            } else {
                if (node.tettree >= 0  &&  node.leg == LegLorR.left) {
                    heights.add(node.height);
                }
            }
        }
    }


    public void moveHybridHeight(int foot1, int foot2, double minh) {
        int f1anc = dhnodes[foot1].anc;
        int f2anc = dhnodes[foot2].anc;
        double maxh = Math.min(dhnodes[f1anc].height, dhnodes[f2anc].height);
        double oldh = dhnodes[foot1].height;
        double newh = AlloppMisc.uniformInRange(oldh, minh, maxh, 0.3);
        dhnodes[foot1].height = dhnodes[foot2].height = newh;
    }


    int scaleAllHeights(double scale) {
        int nofnodes = dhnodes.length;
        int count = 0;
        for (int n = 0; n < nofnodes; n++) {
            if (dhnodes[n].lft >= 0 || dhnodes[n].tettree >= 0)   {
                dhnodes[n].height *= scale;
                count++;
            }
        }
        return count;
    }


    // for assert tests during merging tetrees in move
    public void clearAllNodeTettree() {
        for (int n = 0; n < dhnodes.length; ++n) {
            dhnodes[n].tettree = -1;
        }
    }

    // for moving hyb time, and for assert tests during merging tetrees in move, and
    public int getNodeTettree(int node) {
        return dhnodes[node].tettree;
    }

    LegLorR getNodeLeg(int node) {
        return dhnodes[node].leg;
    }

    // for move that flips all seqs of tet tree and its legs
    void setNodeLeg(int nn, LegLorR leg) {
        dhnodes[nn].leg = leg;
    }


    // for merging tetrees in move
    public void setNodeTettree(int node, int tt) {
        dhnodes[node].tettree = tt;
    }



    /*
      * For testing
      */
    public boolean diphistOK(boolean diprootisroot) {
        int nroots = 0;
        for (int i = 0; i < dhnodes.length; i++) {
            if (dhnodes[i].anc < 0) {
                nroots++;
            }
        }
        if (nroots != 1) {
            return false;
        }
        for (int i = 0; i < dhnodes.length; i++) {
            int nparents = 0;
            for (int j = 0; j < dhnodes.length; j++) {
                if (dhnodes[j].lft == i) { nparents++; }
                if (dhnodes[j].rgt == i) { nparents++; }
            }
            if (dhnodes[i].anc < 0  &&  nparents != 0) {
                return false;
            }
            if (dhnodes[i].anc >= 0  &&  nparents != 1) {
                return false;
            }
        }

        for (int i = 0; i < dhnodes.length; i++) {
            if (dhnodes[i].getNumber() != i) {
                return false;
            }
        }

        for (int i = 0; i < dhnodes.length; i++) {
            if (dhnodes[i].lft >= 0) {
                if (dhnodes[i].rgt < 0) {
                    return false;
                }
                int lft = dhnodes[i].lft;
                int rgt = dhnodes[i].rgt;
                if (dhnodes[lft].anc != i) {
                    return false;
                }
                if (dhnodes[rgt].anc != i) {
                    return false;
                }
                if (dhnodes[i].height <= dhnodes[lft].height) {
                    return false;
                }
                if (dhnodes[i].height <= dhnodes[rgt].height) {
                    return false;
                }
                if (dhnodes[i].tettree >= 0) {
                    return false;
                }
            } else {
                if (dhnodes[i].tettree >= 0) {
                    if (dhnodes[i].height <= 0) {
                        return false;
                    }
                    if (dhnodes[i].leg != LegLorR.left  &&  dhnodes[i].leg != LegLorR.right) {
                        return false;
                    }
                } else {
                    if (dhnodes[i].height != 0) {
                        return false;
                    }
                }
            }
        }
        if (dhnodes[rootn].anc >= 0) {
            return false;
        }
        ArrayList<Integer> feet = collectFeet();
        for (Integer f1 : feet) {
            for (Integer f2 : feet) {
                if (dhnodes[f1].tettree == dhnodes[f2].tettree) {
                    if (dhnodes[f1].height != dhnodes[f2].height) {
                        return false;
                    }
                }
            }
        }
        if (diprootisroot) {
            if (!gotDipTipInSubtree(dhnodes[rootn].lft)  ||  !gotDipTipInSubtree(dhnodes[rootn].rgt)) {
                return false;
            }
        }


        return true;
    }


    private boolean gotDipTipInSubtree(int nn) {
        if (dhnodes[nn].lft < 0) {
            return (dhnodes[nn].tettree < 0);
        } else {
            return gotDipTipInSubtree(dhnodes[nn].lft)  ||  gotDipTipInSubtree(dhnodes[nn].rgt);
        }
    }



    /*
      * **************************************************************************
      *                      PRIVATE methods
      * **************************************************************************
      */





    /*
      * This builds a new tree in dhnodes[] from the one in tmpnodes.
      * The reason for not just copying the array is that the order of
      * nodes in tmpnodes[] is not postorder. (One could live with the
      * nodes in any order, but it would complicate things elsewhere.)
      */
    private void buildSubtreeFromNodes(AlloppSpeciesNetworkModel apspnet,
                                       DipHistNode[] tmpnodes, int n) {
        if (tmpnodes[n].lft < 0) {
            assert tmpnodes[n].rgt < 0;
            dhnodes[nextn] = new DipHistNode(nextn, tmpnodes[n]);
            int tt = dhnodes[nextn].tettree;
            LegLorR leg = dhnodes[nextn].leg;
            AlloppLeggedTree ttree;
            if (tt >= 0) {
                ttree = apspnet.getTetraploidTree(tt);
                if (ttree.getDiphistLftLeg() == n) {
                    if (leg == LegLorR.left) {
                        ttree.setDiphistLftLeg(nextn);
                    } else {
                        assert leg == LegLorR.right;
                        ttree.setDiphistRgtLeg(nextn);
                    }
                } else {
                    assert apspnet.getTetraploidTree(tt).getDiphistRgtLeg() == n;
                    if (leg == LegLorR.left) {
                        ttree.setDiphistLftLeg(nextn);
                    } else {
                        assert leg == LegLorR.right;
                        ttree.setDiphistRgtLeg(nextn);
                    }
                }
            }
            nextn ++;
        } else {
            assert tmpnodes[n].rgt >= 0;
            buildSubtreeFromNodes(apspnet, tmpnodes, tmpnodes[n].lft);
            int lft = nextn - 1;
            buildSubtreeFromNodes(apspnet, tmpnodes, tmpnodes[n].rgt);
            int rgt = nextn - 1;
            dhnodes[nextn] = new DipHistNode(nextn, tmpnodes[n]);
            dhnodes[nextn].lft = lft;
            dhnodes[lft].anc = nextn;
            dhnodes[nextn].rgt = rgt;
            dhnodes[rgt].anc = nextn;
            nextn ++;
        }
    }



    private void removeTip(int leg, DipHistNode[] tmpnodes) {
        int legsibling;
        int leganc = tmpnodes[leg].anc;
        assert leganc >= 0;
        if (tmpnodes[leganc].lft == leg) {
            legsibling = tmpnodes[leganc].rgt;
        } else {
            assert tmpnodes[leganc].rgt == leg;
            legsibling = tmpnodes[leganc].lft;
        }
        int legancanc = tmpnodes[leganc].anc;
        assert legancanc >= 0;
        if (tmpnodes[legancanc].lft == leganc) {
            tmpnodes[legancanc].lft = legsibling;
        } else {
            assert tmpnodes[legancanc].rgt == leganc;
            tmpnodes[legancanc].rgt = legsibling;
        }
    }

    // for testing
    private SimpleTree makesimpletree() {
        SimpleNode[] snodes = new SimpleNode[dhnodes.length];
        for (int n = 0; n < dhnodes.length; n++) {
            snodes[n] = new SimpleNode();
            snodes[n].setTaxon(null); // I use taxon==null to identify joined leg node when removing hybtips
        }
        makesimplesubtree(snodes, 0, dhnodes[rootn]);
        return new SimpleTree(snodes[dhnodes.length-1]);
    }


    // for testing. for makesimpletree()
    private int makesimplesubtree(SimpleNode[] snodes, int nextsn, DipHistNode dhnode) {
        if (dhnode.lft < 0) {
            Taxon tx = new Taxon(dhnode.taxon.getId());
            snodes[nextsn].setTaxon(tx);
            if (dhnode.tettree >= 0) {
                snodes[nextsn].setAttribute("tettree", dhnode.tettree);
                snodes[nextsn].setAttribute("leg", dhnode.leg);
            }
        } else {
            nextsn = makesimplesubtree(snodes, nextsn, dhnodes[dhnode.lft]);
            int subtree0 = nextsn-1;
            nextsn = makesimplesubtree(snodes, nextsn, dhnodes[dhnode.rgt]);
            int subtree1 = nextsn-1;
            snodes[nextsn].addChild(snodes[subtree0]);
            snodes[nextsn].addChild(snodes[subtree1]);
        }
        snodes[nextsn].setHeight(dhnode.height);
        return nextsn+1;
    }


    // for testing
    private void simpletree2dhtesttree(SimpleNode snode) {
        if (snode.getChildCount() == 2) {
            simpletree2dhtesttree(snode.getChild(0));
            int lft = nextn - 1;
            simpletree2dhtesttree(snode.getChild(1));
            int rgt = nextn - 1;
            dhnodes[nextn].lft = lft; dhnodes[lft].anc = nextn;
            dhnodes[nextn].rgt = rgt; dhnodes[rgt].anc = nextn;
        }
        dhnodes[nextn].height = snode.getHeight();
        dhnodes[nextn].taxon = new Taxon(snode.getTaxon().getId());
        dhnodes[nextn].union = apsp.speciesseqEmptyUnion();
        nextn++;
    }









}
