package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.StochasticSimulator;
import dr.evomodel.epidemiology.TwoPathogenModel;
import dr.evomodel.epidemiology.TwoPathogenModelLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

public class TwoPathogenModelLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TWO_PATHOGEN_MODEL_LIKELIHOOD = "twoPathogenModelLikelihood";
    public static final String TREE_MODEL_ONE = "treeModelOne";
    public static final String TREE_MODEL_TWO = "treeModelTwo";
    public static final String STOCHASTIC_SIMULATOR = "stochasticSimulator";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TwoPathogenModel tpm = (TwoPathogenModel) xo.getChild(TwoPathogenModel.class);
        StochasticSimulator stochasticSimulator = (StochasticSimulator) xo.getChild(StochasticSimulator.class);
        TreeModel treeModelOne = (TreeModel) xo.getChild(TREE_MODEL_ONE).getChild(TreeModel.class);
        TreeModel treeModelTwo = (TreeModel) xo.getChild(TREE_MODEL_TWO).getChild(TreeModel.class);

        return new TwoPathogenModelLikelihood(tpm, stochasticSimulator, treeModelOne, treeModelTwo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns stochastic trajectory likelihood for a two-pathogen compartmental model (Shrestha et al., 2011)" +
                "that is used to parameterize birth-death process priors. ";
    }

    public Class getReturnType() {
        return TwoPathogenModelLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserName() {
        return TWO_PATHOGEN_MODEL_LIKELIHOOD;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TwoPathogenModel.class),
            new ElementRule(StochasticSimulator.class),
            new ElementRule(TREE_MODEL_ONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(TreeModel.class),
                    }),
            new ElementRule(TREE_MODEL_TWO,
                    new XMLSyntaxRule[]{
                            new ElementRule(TreeModel.class),
                    }),
    };
}