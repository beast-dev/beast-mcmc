package dr.inferencexml.model;

import dr.inference.model.BooleanLikelihood;
import dr.inference.model.BooleanStatistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

/**
 * Reads a distribution likelihood from a DOM Document element.
 */
public class BooleanLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BOOLEAN_LIKELIHOOD = "booleanLikelihood";

    public static final String DATA = "data";

    public String getParserName() { return BOOLEAN_LIKELIHOOD; }

    public Object parseXMLObject(XMLObject xo) {

        BooleanLikelihood likelihood = new BooleanLikelihood();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof BooleanStatistic) {
                likelihood.addData( (BooleanStatistic)xo.getChild(i));
            }
        }

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A function that log likelihood of a set of boolean statistics. "+
                "If all the statistics are true then it returns 0.0 otherwise -infinity.";
    }

    public Class getReturnType() { return BooleanLikelihood.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(BooleanStatistic.class, 1, Integer.MAX_VALUE )
    };

}
