/*
 * AlloppNetworkNodeSlide.java
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

package dr.evomodel.alloppnet.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.alloppnet.tree.SlidableTree;
import dr.evomodel.alloppnet.speciation.AlloppDiploidHistory;
import dr.evomodel.alloppnet.speciation.AlloppLeggedTree;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.alloppnet.parsers.AlloppNetworkNodeSlideParser;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import jebl.util.FixedBitSet;

import java.util.ArrayList;


/**
 *
 * @author Graham Jones
 *         Date: 01/07/2011
 */


/*
 * An operator for with allopolyploid networks. Uses mnlCanonical() and mnlReconstruct()
 * in SlidableTree.Utils to do some of the work. These are called to change one node in
 * diploid history or one node in a tetratree. There is also a move to change hyb height.
 * 
 * 
 */
public class AlloppNetworkNodeSlide extends SimpleMCMCOperator {

    private final AlloppSpeciesNetworkModel apspnet;
    private final AlloppSpeciesBindings apsp;

    public AlloppNetworkNodeSlide(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double weight) {
        this.apspnet = apspnet;
        this.apsp = apsp;
        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return "None";
    }

    @Override
    public String getOperatorName() {
        return AlloppNetworkNodeSlideParser.NETWORK_NODE_REHEIGHT + "(" + apspnet.getId() +
                "," + apsp.getId() + ")";
    }

    @Override
    public double doOperation() {
        operateOneNodeInNet(0.0);
        return 0;
    }




    private class NodeHeightInNetIndex {
        public int ploidy;  // 2 for diphist, 4 for a tetra tree
        public int tree;    // indexes tettree
        public int index;    // internal node in diphist or tettree, or foot index for hyb height
        public boolean doHybheight;

        public NodeHeightInNetIndex(int ploidy, int tree, int index, boolean doHybheight) {
            this.ploidy = ploidy;
            this.tree = tree;
            this.index = index;
            this.doHybheight = doHybheight;
        }
    }


    private NodeHeightInNetIndex randomnode() {
        int noftettrees = apspnet.getNumberOfTetraTrees();
        int dhcount;
        int hybhcount;
        int count = 0;

        dhcount = apspnet.getNumberOfInternalNodesInDipHist();
        count += dhcount;
        hybhcount = noftettrees;
        count += hybhcount;
        // For each tetratree, the internal/root heights
        for (int i = 0; i < noftettrees; i++) {
            int n = apspnet.getNumberOfInternalNodesInTetTree(i);
            count += n;
        }
        int which = MathUtils.nextInt(count);
        if (which < dhcount) {
            return new NodeHeightInNetIndex(2, 0, which, false);
        } else {
            which -= dhcount;
            if (which < hybhcount) {
                // twice as many feet as hybridizations
                int w = which + (MathUtils.nextBoolean() ? 0 : hybhcount);
                return new NodeHeightInNetIndex(2, 0, w, true);
            } else {
                which -= hybhcount;
                for (int i = 0; i < noftettrees; i++) {
                    int n = apspnet.getNumberOfInternalNodesInTetTree(i);
                    if (which < n) {
                        return new NodeHeightInNetIndex(4, i, which, false);
                    } else {
                        which -= n;
                    }
                }
            }
        }
        assert false;
        return new NodeHeightInNetIndex(-1, -1, -1, false);
    }




    private void operateOneNodeInNet(double factor) {
        assert apspnet.getDiploidHistory().diphistOK(apspnet.getDiploidRootIsRoot());
        NodeHeightInNetIndex nhi = randomnode();
        if (nhi.doHybheight) {
            operateHybridHeight(nhi.index);
        } else {
            if (nhi.ploidy == 2) {
                operateOneNodeInDiploidHistory(nhi.index, factor);
            } else {
                assert nhi.ploidy == 4;
                AlloppLeggedTree altree = apspnet.getTetraploidTree(nhi.tree);
                operateOneNodeInTetraTree(altree, nhi.index, factor);
            }
        }
    }



    private void operateHybridHeight(int footindex) {
        AlloppDiploidHistory diphist = apspnet.getDiploidHistory();
        ArrayList<Integer> feet = diphist.collectFeet();

        assert footindex < feet.size();
        int foot = feet.get(footindex);
        int tt = diphist.getNodeTettree(foot);
        AlloppLeggedTree tettree = apspnet.getTetraploidTree(tt);
        double minh = tettree.getRootHeight();
        int f1 = tettree.getDiphistLftLeg();
        int f2 = tettree.getDiphistRgtLeg();
        assert (foot == f1) || (foot == f2);
        apspnet.beginNetworkEdit();
        diphist.moveHybridHeight(f1, f2, minh);
        apspnet.endNetworkEdit();
    }


    private void operateOneNodeInTetraTree(AlloppLeggedTree tettree, int which, double factor) {

        // As TreeNodeSlide(). Randomly flip children at each node,
        // keeping track of node order (in-order order, left to right).

        NodeRef[] order = SlidableTree.Utils.mnlCanonical(tettree);

        // Find the time of the most recent gene coalescence which
        // has (species,sequence)'s to left and right of this node.
        FixedBitSet left = apsp.speciesseqEmptyUnion();
        FixedBitSet right = apsp.speciesseqEmptyUnion();
        for (int k = 0; k < 2 * which + 1; k += 2) {
            FixedBitSet left0 = apsp.taxonseqToTipUnion(tettree.getSlidableNodeTaxon(order[k]), 0);
            FixedBitSet left1 = apsp.taxonseqToTipUnion(tettree.getSlidableNodeTaxon(order[k]), 1);
            left.union(left0);
            left.union(left1);
        }
        for (int k = 2 * (which + 1); k < order.length; k += 2) {
            FixedBitSet right0 = apsp.taxonseqToTipUnion(tettree.getSlidableNodeTaxon(order[k]), 0);
            FixedBitSet right1 = apsp.taxonseqToTipUnion(tettree.getSlidableNodeTaxon(order[k]), 1);
            right.union(right0);
            right.union(right1);
        }
        double genelimit = apsp.spseqUpperBound(left, right);

        // also keep this node more recent than the hybridization event that led to this tree.
        AlloppDiploidHistory diphist = apspnet.getDiploidHistory();
        double hybridheight = diphist.getHybHeight(tettree);

        final double limit = Math.min(genelimit, hybridheight);

        // On direct call, factor==0.0 and use limit. Else use passed in scaling factor
        double newHeight = -1.0;
        if( factor > 0 ) {
            newHeight = tettree.getSlidableNodeHeight(order[2*which+1]) * factor;
        } else {
            newHeight = MathUtils.nextDouble() * limit;
        }

        apspnet.beginNetworkEdit();
        final NodeRef node = order[2 * which + 1];
        tettree.setSlidableNodeHeight(node, newHeight);
        SlidableTree.Utils.mnlReconstruct(tettree, order);
        apspnet.endNetworkEdit();
    }



    private class RootHeightRange {
        public double lowerlimit;
        public double upperlimit;
        RootHeightRange(double lowerlimit, double upperlimit) {
            this.lowerlimit = lowerlimit;
            this.upperlimit = upperlimit;
        }
    }



    // find limit to keep root a diploid
    // 1. If node to slide is the root, and the second highest node is to left or
    // right of all diploids, then the root must stay the root: lowerlimit =  second highest.
    // 2. If node to slide is not the root, and is to left or right of all diploids,
    // then it must not become the root: upperlimit = root height.
    RootHeightRange findRootRangeForDiploidRootIsRoot(AlloppDiploidHistory diphist, NodeRef[] order, int slidingn) {
        RootHeightRange rootrange = new RootHeightRange(0.0, Double.MAX_VALUE);
        int rootn = -1;
        double maxhgt = 0.0;
        for (int k = 1;  k < order.length;  k += 2) {
            double hgt = diphist.getSlidableNodeHeight(order[k]);
            if (hgt > maxhgt) {
                maxhgt = hgt;
                rootn = k;
            }
        }
        int secondn = -1;
        double secondhgt = 0.0;
        for (int k = 1;  k < order.length;  k += 2) {
            if (k != rootn) {
                double hgt = diphist.getSlidableNodeHeight(order[k]);
                if (hgt > secondhgt) {
                    secondhgt = hgt;
                    secondn = k;
                }
            }
        }
        int leftmostdip = -1;
        int rightmostdip = -1;
        for (int k = 0;  k < order.length;  k += 2) {
            if (diphist.tipIsDiploidTip(order[k])) {
                if (leftmostdip < 0) {
                    leftmostdip = k;
                }
                rightmostdip = k;
            }
        }
        if (slidingn == rootn  &&  (secondn < leftmostdip  ||  secondn > rightmostdip)) {
            rootrange.lowerlimit = diphist.getSlidableNodeHeight(order[secondn]);
        }
        if (slidingn < leftmostdip  ||  slidingn > rightmostdip) {
            rootrange.upperlimit = diphist.getSlidableNodeHeight(order[rootn]);
        }
        return rootrange;

    }


    private void operateOneNodeInDiploidHistory(int which, double factor) {
        apspnet.beginNetworkEdit();
        int slidingn =  2 * which + 1;
        AlloppDiploidHistory diphist = apspnet.getDiploidHistory();

        NodeRef[] order = SlidableTree.Utils.mnlCanonical(diphist);

        // Find the time of the most recent gene coalescence which
        // has (species,sequence)'s to left and right of this node.
        FixedBitSet left = apsp.speciesseqEmptyUnion();
        FixedBitSet right = apsp.speciesseqEmptyUnion();
        for (int k = 0;  k < slidingn;  k += 2) {
            FixedBitSet u = apspnet.calculateDipHistTipUnion(order[k]);
            left.union(u);
        }
        for (int k = slidingn + 1;  k < order.length;  k += 2) {
            FixedBitSet u = apspnet.calculateDipHistTipUnion(order[k]);
            right.union(u);
        }
        double genelimit = apsp.spseqUpperBound(left, right);

        // find limit due to hyb-tips - must be bigger than adjacent heights
        // Note that adjacent nodes in order[] are tips; if they are not hyb-tips
        // they have height zero anyway.
        double hybtiplimit = 0.0;
        if (slidingn-1 >= 0) {
            hybtiplimit = Math.max(hybtiplimit, diphist.getSlidableNodeHeight(order[slidingn-1]));
        }
        if (slidingn+1 < order.length) {
            hybtiplimit = Math.max(hybtiplimit, diphist.getSlidableNodeHeight(order[slidingn+1]));
        }
        RootHeightRange rootrange = new RootHeightRange(0.0, Double.MAX_VALUE);
        if (apspnet.getDiploidRootIsRoot()) {
            rootrange = findRootRangeForDiploidRootIsRoot(diphist, order, slidingn);
        }
        final double upperlimit = Math.min(genelimit, rootrange.upperlimit);
        final double lowerlimit = Math.max(hybtiplimit, rootrange.lowerlimit);

        // On direct call, factor==0.0 and use limit. Else use passed in scaling factor
        double newHeight = -1.0;
        if( factor > 0 ) {
            newHeight = diphist.getSlidableNodeHeight(order[slidingn]) * factor;
        } else {
            newHeight = MathUtils.uniform(lowerlimit, upperlimit);
        }

        assert diphist.diphistOK(apspnet.getDiploidRootIsRoot());
        final NodeRef node = order[slidingn];
        diphist.setSlidableNodeHeight(node, newHeight);
        SlidableTree.Utils.mnlReconstruct(diphist, order);
        if (!diphist.diphistOK(apspnet.getDiploidRootIsRoot())) {
            System.out.println("BUG in operateOneNodeInDiploidHistory()");
        }
        assert diphist.diphistOK(apspnet.getDiploidRootIsRoot());
        apspnet.endNetworkEdit();
    }

}
