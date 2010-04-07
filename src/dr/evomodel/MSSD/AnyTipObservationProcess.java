package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

/**
 * Package: AnyTipObservationProcess
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 18, 2008
 * Time: 6:45:00 PM
 */
public class AnyTipObservationProcess extends AbstractObservationProcess {
    protected double[] u0;
    protected double[] p;

    public AnyTipObservationProcess(String modelName, TreeModel treeModel, PatternList patterns, SiteModel siteModel,
                                    BranchRateModel branchRateModel, Parameter mu, Parameter lam) {
        super(modelName, treeModel, patterns, siteModel, branchRateModel, mu, lam);
    }

    public double calculateLogTreeWeight() {
        int L = treeModel.getNodeCount();
        if (u0 == null || p == null) {
            u0 = new double[L];    // probability that the trait at node i survives to no leaf
            p = new double[L];     // probability of survival on the branch ancestral to i
        }
        int i, j, childNumber;
        NodeRef node;
        double logWeight = 0.0;

        for (i = 0; i < L; ++i) {
            p[i] = 1.0 - getNodeSurvivalProbability(i);
        }

        for (i = 0; i < treeModel.getExternalNodeCount(); ++i) {
            u0[i] = 0.0;
            logWeight += 1.0 - p[i];
        }
        // TODO There is a bug here; the code below assumes nodes are numbered in post-order
        for (i = treeModel.getExternalNodeCount(); i < L; ++i) {
            u0[i] = 1.0;
            node = treeModel.getNode(i);
            for (j = 0; j < treeModel.getChildCount(node); ++j) {
                //childNode = treeModel.getChild(node,j);
                childNumber = treeModel.getChild(node, j).getNumber();
                u0[i] *= 1.0 - p[childNumber] * (1.0 - u0[childNumber]);
            }
            logWeight += (1.0 - u0[i]) * (1.0 - p[i]);
        }

        return -logWeight * lam.getParameterValue(0) / (getAverageRate() * mu.getParameterValue(0));
    }


    private void setTipNodePatternInclusion() { // These values never change
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            final int nChildren = treeModel.getChildCount(node);
            if (nChildren == 0) {
                for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                    extantInTipsBelow[i][patternIndex] = 1;
                    int taxonIndex = patterns.getTaxonIndex(treeModel.getNodeTaxon(node));
                    int[] states = dataType.getStates(patterns.getPatternState(taxonIndex, patternIndex));
                    for (int state : states) {
                        if (state == deathState) {
                            extantInTipsBelow[i][patternIndex] = 0;
                        }
                    }
                    extantInTips[patternIndex] += extantInTipsBelow[i][patternIndex];
                }
            }
        }
    }

    void setNodePatternInclusion() {
        int patternIndex, i, j;
        
        if (nodePatternInclusion == null) {
            nodePatternInclusion = new boolean[nodeCount][patternCount];
        }

        if (this.extantInTips == null) {
            extantInTips = new int[patternCount];
            extantInTipsBelow = new int[nodeCount][patternCount];
            setTipNodePatternInclusion();
        }

        for (patternIndex = 0; patternIndex < patternCount; ++patternIndex) {
            for (i = 0; i < treeModel.getNodeCount(); ++i) {
                NodeRef node = treeModel.getNode(i);
                int nChildren = treeModel.getChildCount(node);
                // TODO There is a bug here; the code below assumes nodes are numbered in post-order
                if (nChildren > 0) {
                    extantInTipsBelow[i][patternIndex] = 0;
                    for (j = 0; j < nChildren; ++j) {
                        int childIndex = treeModel.getChild(node, j).getNumber();
                        extantInTipsBelow[i][patternIndex] += extantInTipsBelow[childIndex][patternIndex];
                    }
                }
            }

            for (i = 0; i < treeModel.getNodeCount(); ++i) {
                nodePatternInclusion[i][patternIndex] = (extantInTipsBelow[i][patternIndex] >= this.extantInTips[patternIndex]);
            }

        }
        nodePatternInclusionKnown = true;
    }

    private int[] extantInTips;
    private int[][] extantInTipsBelow;

}
