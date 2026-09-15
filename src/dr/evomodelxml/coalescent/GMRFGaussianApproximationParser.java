package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFGaussianApproximation;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.xml.*;

public class GMRFGaussianApproximationParser extends AbstractXMLObjectParser {
    public static final String GMRF_GAUSSIAN_APPROXIMATION = "gmrfGaussianApproximation";
    public static final String NEWTON_TOLERANCE = "newtonTolerance";
    public static final String MAX_NEWTON_ITERATIONS = "maxNewtonIterations";

    private static final double DEFAULT_TOLERANCE = 1e-8;
    private static final int DEFAULT_MAX_NEWTON_ITERATIONS = 50;

    public String getParserName(){
        return GMRF_GAUSSIAN_APPROXIMATION;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        GMRFMultilocusSkyrideLikelihood likelihood =
                (GMRFMultilocusSkyrideLikelihood) xo.getChild(GMRFMultilocusSkyrideLikelihood.class);
        if(likelihood == null){
            throw new XMLParseException("gmrfGaussianApproximation requires GMRFMultilocusSkyrideLikelihood element.");
        }

        double tolerance = xo.getAttribute(NEWTON_TOLERANCE, DEFAULT_TOLERANCE);
        int maxIterations = xo.getAttribute(MAX_NEWTON_ITERATIONS, DEFAULT_MAX_NEWTON_ITERATIONS);

        return new GMRFGaussianApproximation(likelihood, tolerance, maxIterations);
    }

    public String getParserDescription(){
        return "Gaussian approximation of skygrid-GLM full conditional distribution of log effective population size vector.";
    }

    public Class getReturnType(){
        return GMRFGaussianApproximation.class;
    }

    public XMLSyntaxRule[] getSyntaxRules(){
        return rules;
    }

    public XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(NEWTON_TOLERANCE, true),
            AttributeRule.newIntegerRule(MAX_NEWTON_ITERATIONS, true),
            new ElementRule(GMRFMultilocusSkyrideLikelihood.class)
    };
}