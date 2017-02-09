package dr.evomodelxml.coalescent.operators;

import dr.evomodel.coalescent.BNPRLikelihood;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.operators.BNPRBlockUpdateOperator;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.IOException;
import java.util.logging.*;


/**
 * Created by mkarcher on 11/3/16.
 */
public class BNPRBlockUpdateOperatorParser extends AbstractXMLObjectParser {
    public static final String BNPR_BLOCK_OPERATOR = "bnprBlockUpdateOperator";

    public static final String KEEP_LOG_RECORD = "keepLogRecord";

    public String getParserName() {
        return BNPR_BLOCK_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean logRecord = xo.getAttribute(KEEP_LOG_RECORD, false);

        Handler bnprHandler;
        Logger bnprLogger = Logger.getLogger("dr.evomodel.coalescent.operators.BNPRBlockUpdateOperator");
        bnprLogger.setUseParentHandlers(false);

        if (logRecord) {
            bnprLogger.setLevel(Level.FINE);

            try {
                bnprHandler = new FileHandler("GMRFBlockUpdate.log." + MathUtils.getSeed());
            } catch (IOException e) {
                throw new RuntimeException(e.getMessage());
            }
            bnprHandler.setLevel(Level.FINE);

            bnprHandler.setFormatter(new XMLFormatter() {
                public String format(LogRecord record) {
                    return "<record>\n \t<message>\n\t" + record.getMessage()
                            + "\n\t</message>\n<record>\n";
                }
            });

            bnprLogger.addHandler(bnprHandler);
        }

//        CoercionMode mode = CoercionMode.parseMode(xo);
//        if (mode == CoercionMode.DEFAULT) mode = CoercionMode.COERCION_ON;

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
//        double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
//        if (scaleFactor == 1.0) {
//            mode = CoercionMode.COERCION_OFF;
//        }

//            if (scaleFactor <= 0.0) {
//                throw new XMLParseException("scaleFactor must be greater than 0.0");
//        if (scaleFactor < 1.0) {
//            throw new XMLParseException("scaleFactor must be greater than or equal to 1.0");
//        }

//        int maxIterations = xo.getAttribute(MAX_ITERATIONS, 200);

//        double stopValue = xo.getAttribute(STOP_VALUE, 0.01);


        BNPRLikelihood bnprLikelihood = (BNPRLikelihood) xo.getChild(BNPRLikelihood.class);
        return new BNPRBlockUpdateOperator(bnprLikelihood, weight);
//        if (xo.getAttribute(OLD_SKYRIDE, true)
//                && !(xo.getName().compareTo(GRID_BLOCK_UPDATE_OPERATOR) == 0)
//                ) {
//            GMRFSkyrideLikelihood gmrfLikelihood = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);
//            return new GMRFSkyrideBlockUpdateOperator(gmrfLikelihood, weight, mode, scaleFactor,
//                    maxIterations, stopValue);
//        } else {
//            GMRFMultilocusSkyrideLikelihood gmrfMultilocusLikelihood = (GMRFMultilocusSkyrideLikelihood) xo.getChild(GMRFMultilocusSkyrideLikelihood.class);
//            return new GMRFMultilocusSkyrideBlockUpdateOperator(gmrfMultilocusLikelihood, weight, mode, scaleFactor,
//                    maxIterations, stopValue);
//        }

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a BNPR block-update operator for the joint distribution of the population sizes and precision parameter.";
    }

    public Class getReturnType() {
        return BNPRBlockUpdateOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
//            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
//            AttributeRule.newDoubleRule(STOP_VALUE, true),
//            AttributeRule.newIntegerRule(MAX_ITERATIONS, true),
//            AttributeRule.newBooleanRule(OLD_SKYRIDE, true),
            new ElementRule(GMRFSkyrideLikelihood.class)
    };
}
