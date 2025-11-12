package dr.inferencexml.distribution;

import dr.inference.distribution.DirichletDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class DirichletDistributionModelParser extends AbstractXMLObjectParser {

    public static final String DIRICHLET_DISTRIBUTION_MODEL = "dirichletDistributionModel";
    public static final String COUNTS = "counts";
    public static final String SUM_TO_NUMBER_OF_ELEMENTS = "sumToNumberOfElements";
    public static final String DISPERSION = "dispersion";

    public String getParserName() {
        return DIRICHLET_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter counts =  (Parameter)xo.getElementFirstChild(COUNTS);
        boolean sumToNumberOfElements = xo.getAttribute(SUM_TO_NUMBER_OF_ELEMENTS, false);

        Parameter dispersion = (Parameter)xo.getElementFirstChild(DISPERSION);

        return new DirichletDistributionModel(counts, dispersion, sumToNumberOfElements);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(COUNTS,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),
            new ElementRule(DISPERSION,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),
            AttributeRule.newBooleanRule(SUM_TO_NUMBER_OF_ELEMENTS, true),
    };

    public String getParserDescription() {
        return "Dirichlet Distribution Model";
    }

    public Class getReturnType() {
        return DirichletDistributionModel.class;
    }


}