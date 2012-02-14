package dr.inferencexml.model;

import dr.xml.*;
import dr.inference.model.Statistic;
import dr.inference.model.OneOnX3Prior;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser of OneOnX3Prior
 *
 */


public class OneOnX3PriorParser extends AbstractXMLObjectParser {

    public static final String ONE_ONE_X_3_PRIOR = "oneOnX3Prior";
    public static final String DATA = "data";

    public String getParserName() {
        return ONE_ONE_X_3_PRIOR;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        OneOnX3Prior likelihood = new OneOnX3Prior();

        XMLObject cxo = xo;

        if (xo.hasChildNamed(DATA)) {
            cxo = xo.getChild(DATA);
        }

        for (int i = 0; i < cxo.getChildCount(); i++) {
            if (cxo.getChild(i) instanceof Statistic) {
                likelihood.addData((Statistic) cxo.getChild(i));
            }
        }

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(Statistic.class, 1, Integer.MAX_VALUE),
                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)})
            )
    };

    public String getParserDescription() {
        return "Calculates the (improper) prior proportional to Prod_i (1/x_i^3) for the given statistic x.";
    }

    public Class getReturnType() {
        return OneOnX3Prior.class;
    }
}

