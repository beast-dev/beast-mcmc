package dr.inferencexml;

import dr.inference.loggers.Logger;
import dr.inference.ml.MLOptimizer;
import dr.inference.model.Likelihood;
import dr.inference.operators.OperatorSchedule;
import dr.xml.*;

import java.util.ArrayList;

/**
 *
 */
public class MLOptimizerParser extends AbstractXMLObjectParser {

    public static final String CHAIN_LENGTH = "chainLength";
    public static final String OPTIMIZER = "optimizer";    

    public String getParserName() { return OPTIMIZER; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int chainLength = xo.getIntegerAttribute(CHAIN_LENGTH);

        OperatorSchedule opsched = null;
        dr.inference.model.Likelihood likelihood = null;
        ArrayList<Logger> loggers = new ArrayList<Logger>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof dr.inference.model.Likelihood) {
                likelihood = (dr.inference.model.Likelihood)child;
            } else if (child instanceof OperatorSchedule) {
                opsched = (OperatorSchedule)child;
            } else if (child instanceof Logger) {
                loggers.add((Logger)child);
            } else {
                throw new XMLParseException("Unrecognized element found in optimizer element:" + child);
            }
        }

        Logger[] loggerArray = new Logger[loggers.size()];
        loggers.toArray(loggerArray);

        return new MLOptimizer("optimizer1", chainLength, likelihood, opsched, loggerArray);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************


    public String getParserDescription() {
        return "This element returns a maximum likelihood heuristic optimizer and runs the optimization as a side effect.";
    }

    public Class getReturnType() { return MLOptimizer.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        AttributeRule.newIntegerRule(CHAIN_LENGTH),
        new ElementRule(OperatorSchedule.class ),
        new ElementRule(Likelihood.class ),
        new ElementRule(Logger.class, 1, Integer.MAX_VALUE )
    };

}
