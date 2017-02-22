
/*
 * AlloppChangeNumHybridizations.java
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

import java.util.ArrayList;


import dr.evomodel.alloppnet.speciation.*;
import dr.evomodel.alloppnet.parsers.AlloppChangeNumHybridizationsParser;
import dr.inference.operators.SimpleMCMCOperator;

import dr.math.MathUtils;
import jebl.util.FixedBitSet;

/**
 * @author Graham Jones
 * Date: 22/07/2012
 */



public class AlloppChangeNumHybridizations  extends SimpleMCMCOperator {

    private final AlloppSpeciesNetworkModel apspnet;
    private final AlloppSpeciesBindings apsp;
    static private final int footdistribution = 2;


    public AlloppChangeNumHybridizations(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double weight) {
        this.apspnet = apspnet;
        this.apsp = apsp;
        setWeight(weight);
    }


    public String getPerformanceSuggestion() {
        return "None";
    }

    @Override
    public String getOperatorName() {
        return AlloppChangeNumHybridizationsParser.CHANGE_NUM_HYBRIDIZATIONS + "(" + apspnet.getId() +
                "," + apsp.getId() + ")";
    }

    @Override
    public double doOperation() {
        if (apspnet.getOneHybridization()) {
            throw new RuntimeException("oneHybridization is true but changeNumHybridizations() called");
        }
        apspnet.beginNetworkEdit();
        double hr = 0.0;
        if (MathUtils.nextBoolean()) {
            hr = doMergeMove();
        } else {
            hr = doSplitMove();
        }
        apspnet.endNetworkEdit();
        assert apspnet.netAndGTreesAreCompatible();
        return hr;
    }

    private class MergeCandidate {
        public int i;
        public int j;

        MergeCandidate(int i, int j) {
            this.i = i;
            this.j = j;
        }
    }


    private class SplitCandidate {
        public int i;
        public AlloppNode root1;
        public AlloppNode root2;

        SplitCandidate(int i, AlloppNode root1, AlloppNode root2) {
            this.i = i;
            this.root1 = root1;
            this.root2 = root2;
        }
    }


    private double doMergeMove() {
        double hr = 0.0;
        ArrayList<MergeCandidate> mcands = findCandidateMerges();
        int nmerges = mcands.size();
        if (nmerges > 0) {
            int mpair = MathUtils.nextInt(nmerges);
            MergeCandidate mcand = mcands.get(mpair);
            hr += Math.log(nmerges);
            hr += mergeTettreePair(mcand.i, mcand.j);
            hr -= Math.log(countCandidateSplits());
            double logpdfoldval = apspnet.removeHybPopParam();
            hr += logpdfoldval;
        }
        return hr;
    }


    private double doSplitMove() {
        double hr = 0.0;
        ArrayList<SplitCandidate> scands = findCandidateSplits();
        int nsplits = scands.size();
        if (nsplits > 0) {
            int stt = MathUtils.nextInt(nsplits);
            SplitCandidate scand = scands.get(stt);
            hr += Math.log(nsplits);
            hr += splitTettree(scand.i, scand.root1, scand.root2);
            hr -= Math.log(countCandidateMerges());
            double logpdfnewval = apspnet.addHybPopParam();
            hr -= logpdfnewval;
        }
        return hr;
    }


    private ArrayList<MergeCandidate> findCandidateMerges() {
        ArrayList<MergeCandidate> mcands = new ArrayList<MergeCandidate>();
        int numttrees = apspnet.getNumberOfTetraTrees();
        for (int i = 0;  i < numttrees;  i++) {
            for (int j = 0;  j < numttrees; j++) {
                if (i != j  &&  pairAreMergeable(i, j)) {
                    mcands.add(new MergeCandidate(i, j));
                }
            }
        }
        return mcands;
    }


    private int countCandidateMerges() {
        return findCandidateMerges().size();
    }



    private ArrayList<SplitCandidate> findCandidateSplits() {
        ArrayList<SplitCandidate> scands = new ArrayList<SplitCandidate>();
        int numttrees = apspnet.getNumberOfTetraTrees();
        for (int i = 0;  i < numttrees;  i++) {
            AlloppLeggedTree ttree = apspnet.getTetraploidTree(i);
            if (ttree.getSlidableNodeCount() > 1) {
                AlloppNode lft = ((AlloppNode)ttree.getSlidableRoot()).getChild(0);
                AlloppNode rgt = ((AlloppNode)ttree.getSlidableRoot()).getChild(1);
                scands.add(new SplitCandidate(i, lft, rgt));
                scands.add(new SplitCandidate(i, rgt, lft));
            }
        }
        return scands;
    }

    private int countCandidateSplits() {
        return findCandidateSplits().size();
    }



    private boolean pairAreMergeable(int tt1, int tt2) {
        boolean mergeable = true;
        AlloppLeggedTree ttree1 = apspnet.getTetraploidTree(tt1);
        AlloppLeggedTree ttree2 = apspnet.getTetraploidTree(tt2);
        AlloppDiploidHistory adhist = apspnet.getDiploidHistory();
        // check legs agree and meet as produced by a split move.
        mergeable = mergeable && adhist.tettreesShareLegs(ttree1, ttree2);
        return mergeable;
    }


    private double mergeTettreePair(int tt1, int tt2) {
        double hr = 0.0;
        AlloppLeggedTree ttree1 = apspnet.getTetraploidTree(tt1);
        AlloppLeggedTree ttree2 = apspnet.getTetraploidTree(tt2);
        AlloppDiploidHistory adhist = apspnet.getDiploidHistory();
        // collect height info
        AlloppDiploidHistory.FootAncHeights lftleg2 =
                adhist.intervalOfFootAncestor(ttree2, AlloppDiploidHistory.LegLorR.left);
        AlloppDiploidHistory.FootAncHeights rgtleg2 =
                adhist.intervalOfFootAncestor(ttree2, AlloppDiploidHistory.LegLorR.right);
        // Choose most recent footanc height as root height of merged tree.
        // Account for loss of the other footanc height.
        // Use gene limit on the lost footanc height for hr calculation.
        // grjtodo-soon test the gene limit calculation somehow
        double rooth;
        if (lftleg2.anchgt < rgtleg2.anchgt) {
            rooth  = lftleg2.anchgt;
            FixedBitSet tt1leg1 = apspnet.unionOfWholeTetTree(tt1, 1);
            FixedBitSet tt2leg1 = apspnet.unionOfWholeTetTree(tt2, 1);
            double genelimit = apsp.spseqUpperBound(tt1leg1, tt2leg1);
            double maxfootanchgt = Math.min(genelimit, rgtleg2.ancanchgt);
            hr += Math.log(uniformpdf(rooth, maxfootanchgt));
        } else {
            rooth  = rgtleg2.anchgt;
            FixedBitSet tt1leg0 = apspnet.unionOfWholeTetTree(tt1, 0);
            FixedBitSet tt2leg0 = apspnet.unionOfWholeTetTree(tt2, 0);
            double genelimit = apsp.spseqUpperBound(tt1leg0, tt2leg0);
            double maxfootanchgt = Math.min(genelimit, lftleg2.ancanchgt);
            hr += Math.log(uniformpdf(rooth, maxfootanchgt));
        }
        // account for loss of two old hybhgts
        hr += Math.log(uniformpdf(ttree1.getRootHeight(), rooth));
        hr += Math.log(uniformpdf(ttree2.getRootHeight(), rooth));
        // merge the trees and replace tt2 with result
        AlloppLeggedTree merged = new AlloppLeggedTree(ttree1, ttree2, rooth);
        apspnet.setTetTree(tt2, merged);
        apspnet.removeTetree(tt1);
        // Fix up the links from diploid history.
        // Get rid of old links first, to enable later assertions
        adhist.clearAllNodeTettree();
        for (int i = 0;  i < apspnet.getNumberOfTetraTrees();  i++) {
            AlloppLeggedTree ttree = apspnet.getTetraploidTree(i);
            int dhlftleg = ttree.getDiphistLftLeg();
            assert adhist.getNodeTettree(dhlftleg) == -1;
            adhist.setNodeTettree(dhlftleg, i);
            int dhrgtleg = ttree.getDiphistRgtLeg();
            assert adhist.getNodeTettree(dhrgtleg) == -1;
            adhist.setNodeTettree(dhrgtleg, i);
        }
        // new hybhgt for merged tree
        double maxhybhgt = Math.min(lftleg2.ancanchgt, rgtleg2.ancanchgt);
        double hybght = MathUtils.uniform(rooth, maxhybhgt);
        adhist.setHybridHeight(merged, hybght);
        hr -= Math.log(uniformpdf(rooth, maxhybhgt));
        adhist.removeFeet(apspnet, ttree1);
        return hr;
    }


    private double uniformpdf(double min, double max) {
        double density = 1.0 / (max-min);
        return density;
    }

    private double splitTettree(int tt, AlloppNode root1, AlloppNode root2) {
        double hr = 0.0;
        // collect info from old TetraTree
        AlloppLeggedTree tetTree = apspnet.getTetraploidTree(tt);
        AlloppDiploidHistory adhist = apspnet.getDiploidHistory();
        double rooth = tetTree.getRootHeight();
        int lftleg = tetTree.getDiphistLftLeg();
        int rgtleg = tetTree.getDiphistRgtLeg();
        double lftanchgt = adhist.getAncHeight(lftleg);
        double rgtanchgt = adhist.getAncHeight(rgtleg);
        // account for the hybhgt that will be lost
        hr += Math.log(uniformpdf(rooth, Math.min(lftanchgt, rgtanchgt)));
        // make two new trees
        AlloppLeggedTree tetTree1 = new AlloppLeggedTree(tetTree, root1);
        AlloppLeggedTree tetTree2 = new AlloppLeggedTree(tetTree, root2);
        // tetree2 gets old one's legs, with new height
        tetTree2.setDiphistLftLeg(lftleg);
        tetTree2.setDiphistRgtLeg(rgtleg);
        double hybhgt2 = MathUtils.uniform(tetTree2.getRootHeight(), rooth);
        hr -= Math.log(uniformpdf(tetTree2.getRootHeight(), rooth));
        adhist.setHybridHeight(tetTree2, hybhgt2);
        // remove old and add new ones to list.
        // tetTree2 replaces tetTree, that is, same index, so dip tips stay consistent
        apspnet.setTetTree(tt, tetTree2);
        int tt2 = tt;
        int tt1 = apspnet.addTetTree(tetTree1);
        // new hybhgt for tree1
        double hybhgt1 = MathUtils.uniform(tetTree1.getRootHeight(), rooth);
        hr -= Math.log(uniformpdf(tetTree1.getRootHeight(), rooth));
        // new hgt for a foot anc (other height is rooth)
        // it is constrained by gene trees and existing node height
        if (MathUtils.nextBoolean()) {
            FixedBitSet tt1leg0 = apspnet.unionOfWholeTetTree(tt1, 0);
            FixedBitSet tt2leg0 = apspnet.unionOfWholeTetTree(tt2, 0);
            double genelimit = apsp.spseqUpperBound(tt1leg0, tt2leg0);
            double maxfootanchgt = Math.min(genelimit, lftanchgt);
            double footanchgt = MathUtils.uniform(rooth, maxfootanchgt);
            hr -= Math.log(uniformpdf(rooth, maxfootanchgt));
            adhist.addTwoDipTips(apspnet, tt1, tt2, footanchgt, rooth, hybhgt1);
        } else {
            FixedBitSet tt1leg1 = apspnet.unionOfWholeTetTree(tt1, 1);
            FixedBitSet tt2leg1 = apspnet.unionOfWholeTetTree(tt2, 1);
            double genelimit = apsp.spseqUpperBound(tt1leg1, tt2leg1);
            double maxfootanchgt = Math.min(genelimit, rgtanchgt);
            double footanchgt = MathUtils.uniform(rooth, maxfootanchgt);
            hr -= Math.log(uniformpdf(rooth, maxfootanchgt));
            adhist.addTwoDipTips(apspnet, tt1, tt2, rooth, footanchgt, hybhgt1);
        }
        // grjtodo-soon The only difference between the two states is the time-order of the nodes.
        // Should topologies or histories be counted?
        // Account for left/right choice. This says histories
        hr += Math.log(2.0);
        return hr;
    }



}




