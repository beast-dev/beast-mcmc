package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

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
    public static final String MODEL_NAME = "anyTipObservationProcess";
    final static String DEATH_RATE = "deathRate";
    final static String IMMIGRATION_RATE = "immigrationRate";
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

    void setNodePatternInclusion() {
        int patternIndex, i, j, extantInTips;
        nodePatternInclusion = new boolean[nodeCount][patternCount];
        for (patternIndex = 0; patternIndex < patternCount; ++patternIndex) {
            extantInTips = 0;

            int states[];
            int extantInTipsBelow[] = new int[treeModel.getNodeCount()];
            for (i = 0; i < treeModel.getNodeCount(); ++i) {
                NodeRef node = treeModel.getNode(i);
                int nChildren = treeModel.getChildCount(node);
                if (nChildren == 0) {    // I'm a tip
                    extantInTipsBelow[i] = 1;
                    int taxonIndex = patterns.getTaxonIndex(treeModel.getNodeTaxon(node));
                    states = dataType.getStates(patterns.getPatternState(taxonIndex, patternIndex));
                    for (int state : states) {
                        if (state == deathState) {
                            extantInTipsBelow[i] = 0;
                        }
                    }
                    extantInTips += extantInTipsBelow[i];
                } else {
                    extantInTipsBelow[i] = 0;
                    for (j = 0; j < nChildren; ++j) {
                        int childIndex = treeModel.getChild(node, j).getNumber();
                        extantInTipsBelow[i] += extantInTipsBelow[childIndex];
                    }
                }
            }

            for (i = 0; i < treeModel.getNodeCount(); ++i) {
                nodePatternInclusion[i][patternIndex] = (extantInTipsBelow[i] >= extantInTips);
            }

        }
        nodePatternInclusionKnown = true;
    }

    /* ***************************************
     *  PARSER IMPLEMENTATION
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter mu = (Parameter) xo.getElementFirstChild(DEATH_RATE);
            Parameter lam = (Parameter) xo.getElementFirstChild(IMMIGRATION_RATE);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            PatternList patterns = (PatternList) xo.getChild(PatternList.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
            Logger.getLogger("dr.evomodel.MSSD").info("Creating AnyTipObservationProcess model. Observed traits are assumed to be extant in at least one tip node. Initial mu = " + mu.getParameterValue(0) + " initial lam = " + lam.getParameterValue(0));

            return new AnyTipObservationProcess(MODEL_NAME, treeModel, patterns, siteModel, branchRateModel, mu, lam);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an instance of the AnyTipObservationProcess for ALSTreeLikelihood calculations";
        }

        public Class getReturnType() {
            return AnyTipObservationProcess.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(PatternList.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class),
                new ElementRule(DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
        };

    };

}
