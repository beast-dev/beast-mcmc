package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.operators.TipStateSwapOperator;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TipStateSwapOperatorParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return TipStateSwapOperator.TIP_STATE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AncestralStateBeagleTreeLikelihood treeLikelihood =
                (AncestralStateBeagleTreeLikelihood) xo.getChild(AncestralStateBeagleTreeLikelihood.class);
        final double weight = xo.getDoubleAttribute("weight");
        return new TipStateSwapOperator(treeLikelihood, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator to swap tip states between two random tips.";
    }

    public Class getReturnType() {
        return TipStateSwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule("weight"),
            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
    };
}
