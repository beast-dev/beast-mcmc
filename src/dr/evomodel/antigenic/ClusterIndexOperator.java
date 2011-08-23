package dr.evomodel.antigenic;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * An integer uniform operator for allocation of items to clusters.
 *
 * @author Andrew Rambaut
 * @version $Id: UniformOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class ClusterIndexOperator extends SimpleMCMCOperator {
    public final static String CLUSTER_INDEX_OPERATOR = "clusterIndexOperator";

    private final int clusterCount;
    private int upper;
    private int lower;

    public ClusterIndexOperator(Variable parameter, int lower, int upper, double weight) { // Bounds.Staircase
        this.parameter = parameter;
        this.clusterCount = upper;
        setWeight(weight);
        this.lower = lower;
        this.upper = upper;

        System.err.println("ClusterIndexOperator with cluster count: " + clusterCount);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return (Parameter) parameter;
    }

    /**
     * @return the Variable this operator acts on.
     */
    public Variable getVariable() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        int index = MathUtils.nextInt(parameter.getSize());

        if (parameter instanceof Parameter) {
            // draw from lower to upper, inclusive
            int newValue = MathUtils.nextInt(upper - lower + 1) + lower; // from 0 to n-1, n must > 0,
            ((Parameter) parameter).setParameterValue(index, newValue);
        } else { // Variable<Integer>, Bounds.Staircase

            int upper = ((Variable<Integer>) parameter).getBounds().getUpperLimit(index);
            int lower = ((Variable<Integer>) parameter).getBounds().getLowerLimit(index);
            int newValue = MathUtils.nextInt(upper - lower + 1) + lower; // from 0 to n-1, n must > 0,
            ((Variable<Integer>) parameter).setValue(index, newValue);

        }

        // Hastings ration of (1/N)^2
        return -2.0 * Math.log(clusterCount);
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "uniformInteger(" + parameter.getId() + ")";
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
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

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String toString() {
        return CLUSTER_INDEX_OPERATOR + "(" + parameter.getId() + ")";
    }

    //PRIVATE STUFF

    private Variable parameter = null;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CLUSTER_INDEX_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            int lower = (int) (double) ((Parameter) parameter).getBounds().getLowerLimit(0);

            int upper = (int) (double) ((Parameter) parameter).getBounds().getUpperLimit(0);

            if (upper == lower || lower == (int) Double.NEGATIVE_INFINITY || upper == (int) Double.POSITIVE_INFINITY) {
                throw new XMLParseException(this.getParserName() + " boundaries not found in parameter "
                        + parameter.getId() + " Use operator lower and upper !");
            }

            return new ClusterIndexOperator((Parameter) parameter, lower, upper, weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that picks new parameter values uniformly at random.";
        }

        public Class getReturnType() {
            return ClusterIndexOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(Parameter.class)
        };
    };
}
