package dr.evomodelxml.branchratemodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.ContinuousTraitBranchRateModel;
import dr.evomodel.branchratemodel.DiscreteTraitBranchRateModel;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class DiscreteTraitBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String RATES = "rates";
    public static final String INDICATORS = "indicators";
    public static final String TRAIT_INDEX = "traitIndex";
    public static final String TRAIT_NAME = "traitName";

    public String getParserName() {
        return DiscreteTraitBranchRateModel.DISCRETE_TRAIT_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);

        TreeTraitProvider traitProvider = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);

        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);
        Parameter indicatorsParameter = null;

        if (xo.getChild(INDICATORS) != null) {
            indicatorsParameter = (Parameter) xo.getElementFirstChild(INDICATORS);
        }

        int traitIndex = xo.getAttribute(TRAIT_INDEX, 1) - 1;
        String traitName = xo.getAttribute(TRAIT_NAME, "states");

        Logger.getLogger("dr.evomodel").info("Using discrete trait branch rate model.\n" +
                "\tIf you use this model, please cite:\n" +
                "\t\tDrummond and Suchard (in preparation)");

        if (traitProvider == null) {
            // Use the version that reconstructs the trait using parsimony:
            return new DiscreteTraitBranchRateModel(treeModel, patternList, traitIndex, ratesParameter);
        } else {
            TreeTrait trait = traitProvider.getTreeTrait(traitName);
            if (trait == null) {
                throw new XMLParseException("A trait called, " + traitName + ", was not available from the TreeTraitProvider supplied to " + getParserName() + ", with ID " + xo.getId());
            }

            return new DiscreteTraitBranchRateModel(treeModel, trait, traitIndex, ratesParameter, indicatorsParameter);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This Branch Rate Model takes a discrete trait reconstruction (provided by a TreeTraitProvider) and " +
                        "gives the rate for each branch of the tree based on the child trait of " +
                        "that branch. The rates for each trait value are specified in a multidimensional parameter.";
    }

    public Class getReturnType() {
        return DiscreteTraitBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class, "The tree model"),
            new XORRule(
                    new ElementRule(TreeTraitProvider.class, "The trait provider"),
                    new ElementRule(PatternList.class)),
            new ElementRule(RATES, Parameter.class, "The rates of the different trait values", false),
            new ElementRule(INDICATORS, Parameter.class, "An index that links the state with a rate", true),
            AttributeRule.newIntegerRule(TRAIT_INDEX, true),
            AttributeRule.newStringRule(TRAIT_NAME, true)
    };

}