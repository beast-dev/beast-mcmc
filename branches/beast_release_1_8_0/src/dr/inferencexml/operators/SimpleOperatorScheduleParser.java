package dr.inferencexml.operators;

import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.xml.*;

import java.util.logging.Logger;

/**
 *
 */
public class SimpleOperatorScheduleParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_SCHEDULE = "operators";
    public static final String SEQUENTIAL = "sequential";
    public static final String OPTIMIZATION_SCHEDULE = "optimizationSchedule";

    public String getParserName() {
        return OPERATOR_SCHEDULE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SimpleOperatorSchedule schedule = new SimpleOperatorSchedule();

        if (xo.hasAttribute(SEQUENTIAL)) {
            schedule.setSequential(xo.getBooleanAttribute(SEQUENTIAL));
        }


        if (xo.hasAttribute(OPTIMIZATION_SCHEDULE)) {
            String type = xo.getStringAttribute(OPTIMIZATION_SCHEDULE);
            Logger.getLogger("dr.inference").info("Optimization Schedule: " + type);

            if (type.equals(OperatorSchedule.LOG_STRING))
                schedule.setOptimizationSchedule(OperatorSchedule.LOG_SCHEDULE);
            else if (type.equals(OperatorSchedule.SQRT_STRING))
                schedule.setOptimizationSchedule(OperatorSchedule.SQRT_SCHEDULE);
            else if (!type.equals(OperatorSchedule.DEFAULT_STRING))
                throw new RuntimeException("Unsupported optimization schedule");
        }

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof MCMCOperator) {
                schedule.addOperator((MCMCOperator) child);
            }
        }
        return schedule;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SEQUENTIAL, true),
            new ElementRule(MCMCOperator.class, 1, Integer.MAX_VALUE),
            AttributeRule.newStringRule(OPTIMIZATION_SCHEDULE, true)
    };

    public String getParserDescription() {
        return "A simple operator scheduler";
    }

    public Class getReturnType() {
        return SimpleOperatorSchedule.class;
    }
    
}
