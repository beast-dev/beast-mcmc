package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.PositiveDefiniteSubstitutionModel;
import dr.inference.model.MatrixParameter;
import dr.xml.*;

/**
 *
 */
public class PositiveDefiniteSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String SVS_GENERAL_SUBSTITUTION_MODEL = "positiveDefiniteSubstitutionModel";

    public String getParserName() {
        return SVS_GENERAL_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameter ratesParameter = (MatrixParameter) xo.getChild(MatrixParameter.class);

        return new PositiveDefiniteSubstitutionModel(ratesParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of sequence substitution for any data type with stochastic variable selection.";
    }

    public Class getReturnType() {
        return PositiveDefiniteSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(MatrixParameter.class)
    };

}
