package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.GaussianProcessSkytrackLikelihood;
import dr.evomodel.coalescent.operators.GaussianProcessSkytrackBlockUpdateOperator;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.IOException;
import java.util.logging.*;

/**
 *
 */
public class GaussianProcessSkytrackBlockUpdateOperatorParser extends AbstractXMLObjectParser {

    public static final String BLOCK_UPDATE_OPERATOR = "gpBlockUpdateOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String MAX_ITERATIONS = "maxIterations";
    public static final String STOP_VALUE = "stopValue";
    public static final String KEEP_LOG_RECORD = "keepLogRecord";
    public static final String OLD_SKYRIDE = "oldSkyride";

    public String getParserName() {
        return BLOCK_UPDATE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean logRecord = xo.getAttribute(KEEP_LOG_RECORD, false);

        Handler gmrfHandler;
        Logger gmrfLogger = Logger.getLogger("dr.evomodel.coalescent.operators.GaussianProcessSkytrackBlockUpdateOperator");
        gmrfLogger.setUseParentHandlers(false);

        if (logRecord) {
            gmrfLogger.setLevel(Level.FINE);

            try {
                gmrfHandler = new FileHandler("GPBlockUpdate.log." + MathUtils.getSeed());
            } catch (IOException e) {
                throw new RuntimeException(e.getMessage());
            }
            gmrfHandler.setLevel(Level.FINE);

            gmrfHandler.setFormatter(new XMLFormatter() {
                public String format(LogRecord record) {
                    return "<record>\n \t<message>\n\t" + record.getMessage()
                            + "\n\t</message>\n<record>\n";
                }
            });

            gmrfLogger.addHandler(gmrfHandler);
        }

        CoercionMode mode = CoercionMode.parseMode(xo);
        if (mode == CoercionMode.DEFAULT) mode = CoercionMode.COERCION_ON;

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

//            if (scaleFactor <= 0.0) {
//                throw new XMLParseException("scaleFactor must be greater than 0.0");
        if (scaleFactor < 1.0) {
            throw new XMLParseException("scaleFactor must be greater than or equal to 1.0");
        }

        int maxIterations = xo.getAttribute(MAX_ITERATIONS, 200);

        double stopValue = xo.getAttribute(STOP_VALUE, 0.01);


        if (xo.getAttribute(OLD_SKYRIDE, true)) {

            GaussianProcessSkytrackLikelihood gmrfLikelihood = (GaussianProcessSkytrackLikelihood) xo.getChild(GaussianProcessSkytrackLikelihood.class);
            return new GaussianProcessSkytrackBlockUpdateOperator(gmrfLikelihood, weight, mode, scaleFactor,
                maxIterations, stopValue);
        }else{
           GaussianProcessSkytrackLikelihood gmrfLikelihood = (GaussianProcessSkytrackLikelihood) xo.getChild(GaussianProcessSkytrackLikelihood.class);
            return new GaussianProcessSkytrackBlockUpdateOperator(gmrfLikelihood, weight, mode, scaleFactor,
                maxIterations, stopValue);
        }

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************


    public String getParserDescription() {
        return "This element returns a GMRF block-update operator for the joint distribution of the population sizes and precision parameter.";
    }

    public Class getReturnType() {
         return GaussianProcessSkytrackBlockUpdateOperator.class;
//        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newDoubleRule(STOP_VALUE, true),
            AttributeRule.newIntegerRule(MAX_ITERATIONS, true),
            AttributeRule.newBooleanRule(OLD_SKYRIDE, true),
    new ElementRule(GaussianProcessSkytrackLikelihood.class)
    };

}
