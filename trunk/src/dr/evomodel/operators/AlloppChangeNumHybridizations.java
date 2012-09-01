
package dr.evomodel.operators;

import java.util.ArrayList;


import dr.evomodel.speciation.*;
import dr.evomodelxml.operators.AlloppChangeNumHybridizationsParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;

import dr.math.MathUtils;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 22/07/12
 */
public class AlloppChangeNumHybridizations  extends SimpleMCMCOperator {

    private final AlloppSpeciesNetworkModel apspnet;
    private final AlloppSpeciesBindings apsp;

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
    public double doOperation() throws OperatorFailedException {
        apspnet.beginNetworkEdit();
        double hr = 0.0;
        if (MathUtils.nextBoolean()) {
            hr = doMergeMove();
        } else {
            hr = doSplitMove();
        }
        apspnet.endNetworkEdit();
        assert apspnet.getDiploidHistory().diphistOK();
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
        double hybh1 = adhist.getHybHeight(ttree1);
        double hybh2 = adhist.getHybHeight(ttree2);
        double rooth2 = ttree2.getRootHeight();
        mergeable = mergeable && (rooth2 <= hybh1);
        mergeable = mergeable && (hybh1 <= hybh2);
        mergeable = mergeable && adhist.tettreesShareLegs(ttree1, ttree2);
        return mergeable;
    }


    private double mergeTettreePair(int tt1, int tt2) {
        double hr = 0.0;
        AlloppLeggedTree ttree1 = apspnet.getTetraploidTree(tt1);
        AlloppLeggedTree ttree2 = apspnet.getTetraploidTree(tt2);
        AlloppDiploidHistory adhist = apspnet.getDiploidHistory();

        double lfttime = adhist.intervalOfFoot(ttree2, true);
        double rgttime = adhist.intervalOfFoot(ttree2, false);
        hr -= Math.log(lfttime * rgttime);

        // merge the trees and replace tt2 with result
        AlloppLeggedTree merged = new AlloppLeggedTree(ttree1, ttree2, adhist.getHybHeight(ttree1));
        apspnet.setTetTree(tt2, merged);
        apspnet.removeTetree(tt1);

        // fix up the links from diploid history
        // get rid of old links first, to enable later assertions
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
        adhist.removeFeet(apspnet, ttree1);
        return hr;
    }



    private double splitTettree(int tt, AlloppNode root1, AlloppNode root2) {
        double hr = 0.0;
        AlloppLeggedTree tetTree = apspnet.getTetraploidTree(tt);
        AlloppDiploidHistory adhist = apspnet.getDiploidHistory();
        double rooth = tetTree.getRootHeight();
        int lftleg = tetTree.getDiphistLftLeg();
        int rgtleg = tetTree.getDiphistRgtLeg();
        double hybh = adhist.getHybHeight(tetTree);
        double lftanchgt = adhist.getAncHeight(lftleg);
        double rgtanchgt = adhist.getAncHeight(rgtleg);

        // make two new trees, number 2 getting old one's legs.
        AlloppLeggedTree tetTree1 = new AlloppLeggedTree(tetTree, root1);
        AlloppLeggedTree tetTree2 = new AlloppLeggedTree(tetTree, root2);
        tetTree2.setDiphistLftLeg(tetTree.getDiphistLftLeg());
        tetTree2.setDiphistRgtLeg(tetTree.getDiphistRgtLeg());
        // remove old and add new ones to list.
        // tetTree2 replaces tetTree, that is, same index, so dip tips stay consistent
        apspnet.setTetTree(tt, tetTree2);
        int tt2 = tt;
        int tt1 = apspnet.addTetTree(tetTree1);

        double lfthgt = MathUtils.uniform(hybh, lftanchgt);
        double rgthgt = MathUtils.uniform(hybh, rgtanchgt);
        hr += Math.log((lftanchgt-hybh) * (rgtanchgt-hybh));
        adhist.addTwoDipTips(apspnet, tt1, tt2, lfthgt, rgthgt, rooth);
        return hr;
    }



}




