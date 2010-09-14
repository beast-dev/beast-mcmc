package dr.evomodel.operators;

import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.IntegratedMultivariateTraitLikelihood;
import dr.evomodel.treelikelihood.AbstractLikelihoodCore;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TipTraitSwapOperator extends SimpleMCMCOperator {

    public static final String TIP_TRAIT_SWAP_OPERATOR = "tipTraitSwapOperator";

    public TipTraitSwapOperator(String traitName, IntegratedMultivariateTraitLikelihood traitLikelihood, double weight) {
        this.traitLikelihood = traitLikelihood;
        this.traitName = traitName;
        setWeight(weight);
    }

    private int index1;
    private int index2;

    public double doOperation() throws OperatorFailedException {

        int tipCount = traitLikelihood.getTreeModel().getExternalNodeCount();

        // Choose two tips to swap
        index1 = MathUtils.nextInt(tipCount);
        index2 = index1;
        while (index2 == index1)
            index2 = MathUtils.nextInt(tipCount);

        swap(index1, index2);

        traitLikelihood.makeDirty();

        return 0;
    }

    private void swap(int i, int j) {
        double[] trait1 = traitLikelihood.getTipDataValues(i);
        double[] trait2 = traitLikelihood.getTipDataValues(j);

        traitLikelihood.setTipDataValuesForNode(j, trait1);
        traitLikelihood.setTipDataValuesForNode(i, trait2);           
    }

    public void reject() {
        super.reject();
        // There is currently no restore functions for tip states, so manually adjust state
        swap(index1, index2);
    }

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String getOperatorName() {
        return TIP_TRAIT_SWAP_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TIP_TRAIT_SWAP_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            IntegratedMultivariateTraitLikelihood traitLikelihood =
                    (IntegratedMultivariateTraitLikelihood) xo.getChild(IntegratedMultivariateTraitLikelihood.class);
            final String traitName = xo.getStringAttribute("trait");
            final double weight = xo.getDoubleAttribute("weight");
            return new TipTraitSwapOperator(traitName, traitLikelihood, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an operator to swap tip traits between two random tips.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("trait"),
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(IntegratedMultivariateTraitLikelihood.class),
        };

    };


    private final IntegratedMultivariateTraitLikelihood traitLikelihood;
    private final String traitName;
}
