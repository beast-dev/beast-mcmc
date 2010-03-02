package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFTestLikelihood;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 *
 */
public class GMRFTestLikelihoodParser extends GMRFSkyrideLikelihoodParser {

    public static final String SKYLINE_TEST_LIKELIHOOD = "gmrfTestLikelihood";
    public static final String INTERVAL_PARAMETER = "intervals";
    public static final String SUFFSTAT_PARAMETER = "sufficientStatistics";

    public String getParserName() {
        return SKYLINE_TEST_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(POPULATION_PARAMETER);
        Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(INTERVAL_PARAMETER);
        Parameter intervalParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(SUFFSTAT_PARAMETER);
        Parameter statParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(LAMBDA_PARAMETER);
        Parameter lambda = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(BETA_PARAMETER);
        Parameter betaParameter = (Parameter) cxo.getChild(Parameter.class);

        DesignMatrix designMatrix = (DesignMatrix) xo.getChild(DesignMatrix.class);


        return new GMRFTestLikelihood(popParameter, precParameter, lambda, betaParameter, designMatrix, intervalParameter, statParameter);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return GMRFTestLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),

            new ElementRule(LAMBDA_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
    };
    
}
