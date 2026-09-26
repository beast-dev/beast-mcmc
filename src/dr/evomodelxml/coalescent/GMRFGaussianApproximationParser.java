package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFGaussianApproximation;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.xml.*;

public class GMRFGaussianApproximationParser extends AbstractXMLObjectParser {
    public static final String GMRF_GAUSSIAN_APPROXIMATION = "gmrfGaussianApproximation";
    public static final String NEWTON_TOLERANCE = "newtonTolerance";
    public static final String MAX_NEWTON_ITERATIONS = "maxNewtonIterations";
    public static final String NUM_IMPORTANCE_SAMPLES = "numImportanceSamples";
    public static final String NUM_QUADRATURE_POINTS = "numQuadraturePoints";
    public static final String LOG_TAU_LOWER = "logTauLower";
    public static final String LOG_TAU_UPPER = "logTauUpper";
    public static final String NUM_IMPORTANCE_SAMPLES_PER_NODE = "numImportanceSamplesPerNode";

    private static final double DEFAULT_TOLERANCE = 1e-8;
    private static final int DEFAULT_MAX_NEWTON_ITERATIONS = 50;
    private static final int DEFAULT_NUM_IMPORTANCE_SAMPLES = 200;
    private static final int DEFAULT_NUM_QUADRATURE_POINTS = 401;
    private static final double DEFAULT_LOG_TAU_LOWER = -15.0;
    private static final double DEFAULT_LOG_TAU_UPPER = 15.0;
    private static final int DEFAULT_NUM_IMPORTANCE_SAMPLES_PER_NODE = 200;

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
        int numImportanceSamples = xo.getAttribute(NUM_IMPORTANCE_SAMPLES, DEFAULT_NUM_IMPORTANCE_SAMPLES);
        int numQuadraturePoints = xo.getAttribute(NUM_QUADRATURE_POINTS, DEFAULT_NUM_QUADRATURE_POINTS);
        double logTauLower = xo.getAttribute(LOG_TAU_LOWER, DEFAULT_LOG_TAU_LOWER);
        double logTauUpper = xo.getAttribute(LOG_TAU_UPPER, DEFAULT_LOG_TAU_UPPER);
        int numImportanceSamplesPerNode = xo.getAttribute(NUM_IMPORTANCE_SAMPLES_PER_NODE,
                DEFAULT_NUM_IMPORTANCE_SAMPLES_PER_NODE);

        return new GMRFGaussianApproximation(likelihood, tolerance, maxIterations, numImportanceSamples,
                numQuadraturePoints, logTauLower, logTauUpper, numImportanceSamplesPerNode);
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
            AttributeRule.newIntegerRule(NUM_IMPORTANCE_SAMPLES, true),
            AttributeRule.newIntegerRule(NUM_QUADRATURE_POINTS, true),
            AttributeRule.newDoubleRule(LOG_TAU_LOWER, true),
            AttributeRule.newDoubleRule(LOG_TAU_UPPER, true),
            AttributeRule.newIntegerRule(NUM_IMPORTANCE_SAMPLES_PER_NODE, true),
            new ElementRule(GMRFMultilocusSkyrideLikelihood.class)
    };
}