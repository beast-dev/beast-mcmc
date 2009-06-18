package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic random walk operator for use with a multi-dimensional Integer parameters.
 *
 * @author Michael Defoin Platel
 * @version $Id: RandomWalkIntegerOperator.java$
 */
public class RandomWalkIntegerOperator extends SimpleMCMCOperator {

    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";

    public RandomWalkIntegerOperator(Parameter parameter, int windowSize, double weight) {
        this.parameter = parameter;
        this.windowSize = windowSize;
        setWeight(weight);
    }


    public RandomWalkIntegerOperator(Parameter parameter, Parameter updateIndex, int windowSize, double weight) {
        this.parameter = parameter;
        this.windowSize = windowSize;
        setWeight(weight);

        updateMap = new ArrayList<Integer>();
        for (int i = 0; i < updateIndex.getDimension(); i++) {
            if (updateIndex.getParameterValue(i) == 1.0)
                updateMap.add(i);
        }
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    public final int getWindowSize() {
        return windowSize;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        // a random dimension to perturb
        int index;
        if (updateMap == null)
            index = MathUtils.nextInt(parameter.getDimension());
        else
            index = updateMap.get(MathUtils.nextInt(updateMap.size()));

        // a random non zero integer around old value within windowSize * 2
        int oldValue = (int) parameter.getParameterValue(index);
        int newValue;
        int roll = MathUtils.nextInt(2 * windowSize);
        if (roll >= windowSize) {
            newValue = oldValue + 1 + roll - windowSize;

            if (newValue > parameter.getBounds().getUpperLimit(index))
                newValue = 2 * (int) parameter.getBounds().getUpperLimit(index) - newValue;
        } else {
            newValue = oldValue - 1 - roll;

            if (newValue < parameter.getBounds().getLowerLimit(index))
                newValue = 2 * (int) parameter.getBounds().getLowerLimit(index) - newValue;
        }

        parameter.setParameterValue(index, newValue);

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "randomWalkInteger(" + parameter.getParameterName() + ")";
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(windowSize, parameter.getParameterValue(0) * 2.0, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }

    public static dr.xml.XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return "randomWalkIntegerOperator";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            double d = xo.getDoubleAttribute(WINDOW_SIZE);
            if (d != Math.floor(d)) {
                throw new XMLParseException("The window size of a randomWalkIntegerOperator should be an integer");
            }

            int windowSize = (int)d;
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            if (xo.hasChildNamed(UPDATE_INDEX)) {
                XMLObject cxo = (XMLObject) xo.getChild(UPDATE_INDEX);
                Parameter updateIndex = (Parameter) cxo.getChild(Parameter.class);
                return new RandomWalkIntegerOperator(parameter, updateIndex, windowSize, weight);
            }

            return new RandomWalkIntegerOperator(parameter, windowSize, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WINDOW_SIZE),
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(Parameter.class)
        };

    };

    public String toString() {
        return "randomWalkIntegerOperator(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private int windowSize = 1;
    private List<Integer> updateMap = null;
}
