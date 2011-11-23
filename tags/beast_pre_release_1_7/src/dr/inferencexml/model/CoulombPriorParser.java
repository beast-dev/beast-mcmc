package dr.inferencexml.model;

import dr.inference.model.CoulombPrior;
import dr.inference.model.OneOnXPrior;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * Reads a distribution likelihood from a DOM Document element.
 */
public class CoulombPriorParser extends AbstractXMLObjectParser {

    public static final String COULOMB_PRIOR = "coulombPrior";
    public static final String BETA = "beta";
    public static final String DATA = "data";

    public String getParserName() {
        return COULOMB_PRIOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo;
        final double beta = xo.getDoubleAttribute(BETA);

        if (xo.hasChildNamed(DATA)) {
            cxo = xo.getChild(DATA);
        }

        CoulombPrior likelihood = new CoulombPrior(beta);

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
            AttributeRule.newDoubleRule(BETA),
            new XORRule(
                    new ElementRule(Statistic.class, 1, Integer.MAX_VALUE),
                    new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)})
            )
    };

    public String getParserDescription() {
        return "Calculates a prior density based on the force due to electrostatic charge between particles the given distance x.";
    }

    public Class getReturnType() {
        return OneOnXPrior.class;
    }
}
