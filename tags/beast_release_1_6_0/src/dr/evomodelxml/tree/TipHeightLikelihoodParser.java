package dr.evomodelxml.tree;

import dr.evomodel.tree.TipHeightLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a distribution likelihood from a DOM Document element.
 */
public class TipHeightLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TIP_HEIGHT_LIKELIHOOD = "tipHeightLikelihood";

    public static final String DISTRIBUTION = "distribution";
    public static final String TIP_HEIGHTS = "tipHeights";

    public String getParserName() {
        return TIP_HEIGHT_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ParametricDistributionModel model = (ParametricDistributionModel) xo.getElementFirstChild(DISTRIBUTION);
        Parameter tipHeights = (Parameter) xo.getElementFirstChild(TIP_HEIGHTS);

        return new TipHeightLikelihood(model, tipHeights);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(DISTRIBUTION,
                    new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)}),
            new ElementRule(TIP_HEIGHTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    public String getParserDescription() {
        return "Calculates the likelihood of the tipHeights given some parametric or empirical distribution.";
    }

    public Class getReturnType() {
        return TipHeightLikelihood.class;
    }
}
