/*
 * UpDownOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.operators;

import dr.math.MathUtils;

public class UpDownOperator extends AbstractCoercableOperator {

    private Scalable[] upParameter = null;
    private Scalable[] downParameter = null;
    private double scaleFactor;

    public UpDownOperator(Scalable[] upParameter, Scalable[] downParameter,
                          double scale, double weight, CoercionMode mode) {

        super(mode);
        setWeight(weight);

        this.upParameter = upParameter;
        this.downParameter = downParameter;
        this.scaleFactor = scale;
    }

    public final double getScaleFactor() {
        return scaleFactor;
    }

    public final void setScaleFactor(double sf) {
        if( (sf > 0.0) && (sf < 1.0) ) {
            scaleFactor = sf;
        } else {
            throw new IllegalArgumentException("scale must be between 0 and 1");
        }
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
        int goingUp = 0, goingDown = 0;

        if( upParameter != null ) {
            for( Scalable up : upParameter ) {
                goingUp += up.scale(scale, -1, false);
//                try {
//                    goingUp += up.scale(scale, -1, true);
//                } catch (RuntimeException re) {
//                    return Double.NEGATIVE_INFINITY;
//                }
            }
            for( Scalable up : upParameter ) {
                if (!up.testBounds()) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        if( downParameter != null ) {
            for( Scalable dn : downParameter ) {
                goingDown += dn.scale(1.0 / scale, -1, false);
//                try {
//                    goingDown += dn.scale(1.0 / scale, -1, true);
//                } catch (RuntimeException re) {
//                    return Double.NEGATIVE_INFINITY;
//                }
            }
            for( Scalable dn : downParameter ) {
                if (!dn.testBounds()) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }


        return (goingUp - goingDown - 2) * Math.log(scale);
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public final String getOperatorName() {
        String name = "";
        if( upParameter != null ) {
            name = "up:";
            for( Scalable up : upParameter ) {
                name = name + up.getName() + " ";
            }
        }

        if( downParameter != null ) {
            name += "down:";
            for( Scalable dn : downParameter ) {
                name = name + dn.getName() + " ";
            }
        }
        return name;
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0) / Math.log(10);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.pow(10.0, value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    // Since this operator invariably modifies at least 2 parameters it
    // should allow lower acceptance probabilities
    // as it is known that optimal acceptance levels are inversely
    // proportional to the number of dimensions operated on
    // AD 16/3/2004
    public double getMinimumAcceptanceLevel() {
        return 0.05;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.3;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.10;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.20;
    }
}

// The old implementation for reference and until the new one is tested :)

//
//public class UpDownOperator extends AbstractCoercableOperator {
//
//    public static final String UP_DOWN_OPERATOR = "upDownOperator";
//    public static final String UP = "up";
//    public static final String DOWN = "down";
//
//    public static final String SCALE_FACTOR = "scaleFactor";
//
//    public UpDownOperator(Parameter upParameter, Parameter downParameter,
//                          double scale, double weight, CoercionMode mode) {
//
//        super(mode);
//
//        this.upParameter = upParameter;
//        this.downParameter = downParameter;
//        this.scaleFactor = scale;
//        setWeight(weight);
//    }
//
//    public final int getPriorType() {
//        return prior.getPriorType();
//    }
//
//    public final double getScaleFactor() {
//        return scaleFactor;
//    }
//
//    public final void setPriorType(int type) {
//        prior.setPriorType(type);
//    }
//
//    public final void setScaleFactor(double sf) {
//        if ((sf > 0.0) && (sf < 1.0)) scaleFactor = sf;
//        else throw new IllegalArgumentException("minimum scale must be between 0 and 1");
//    }
//
//    /**
//     * change the parameter and return the hastings ratio.
//     */
//    public final double doOperation() throws OperatorFailedException {
//
//        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
//
//        if( upParameter != null ) {
//            for(int i = 0; i < upParameter.getDimension(); i++) {
//                upParameter.setParameterValue(i, upParameter.getParameterValue(i) * scale);
//            }
//            for(int i = 0; i < upParameter.getDimension(); i++) {
//                if( upParameter.getParameterValue(i) < upParameter.getBounds().getLowerLimit(i) ||
//                        upParameter.getParameterValue(i) > upParameter.getBounds().getUpperLimit(i) ) {
//                    throw new OperatorFailedException("proposed value outside boundaries");
//                }
//            }
//        }
//
//        if( downParameter != null ) {
//            for(int i = 0; i < downParameter.getDimension(); i++) {
//                downParameter.setParameterValue(i, downParameter.getParameterValue(i) / scale);
//            }
//            for(int i = 0; i < downParameter.getDimension(); i++) {
//                if( downParameter.getParameterValue(i) < downParameter.getBounds().getLowerLimit(i) ||
//                        downParameter.getParameterValue(i) > downParameter.getBounds().getUpperLimit(i) ) {
//                    throw new OperatorFailedException("proposed value outside boundaries");
//                }
//            }
//        }
//
//        final int goingUp = (upParameter == null) ? 0 : upParameter.getDimension();
//        final int goingDown = (downParameter == null) ? 0 : downParameter.getDimension();
//
//        final double logq = (goingUp - goingDown - 2) * Math.log(scale);
//
//        return logq;
//    }
//
//    public final String getPerformanceSuggestion() {
//
//        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//        double targetProb = getTargetAcceptanceProbability();
//        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
//        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
//        if (prob < getMinimumGoodAcceptanceLevel()) {
//            return "Try setting scaleFactor to about " + formatter.format(sf);
//        } else if (prob > getMaximumGoodAcceptanceLevel()) {
//            return "Try setting scaleFactor to about " + formatter.format(sf);
//        } else return "";
//    }
//
//    //MCMCOperator INTERFACE
//    public final String getOperatorName() {
//        return (upParameter != null ? "up:" + upParameter.getParameterName() : "") +
//                (downParameter != null ? " down:" + downParameter.getParameterName() : "");
//    }
//
//    public double getCoercableParameter() {
//        return Math.log(1.0 / scaleFactor - 1.0) / Math.log(10);
//    }
//
//    public void setCoercableParameter(double value) {
//        scaleFactor = 1.0 / (Math.pow(10.0, value) + 1.0);
//    }
//
//    public double getRawParameter() {
//        return scaleFactor;
//    }
//
//    public double getTargetAcceptanceProbability() {
//        return 0.234;
//    }
//
//    // Since this operator invariably modifies at least 2 parameters it
//    // should allow lower acceptance probabilities
//    // as it is known that optimal acceptance levels are inversely
//    // proportional to the number of dimensions operated on
//    // AD 16/3/2004
//    public double getMinimumAcceptanceLevel() {
//        return 0.05;
//    }
//
//    public double getMaximumAcceptanceLevel() {
//        return 0.3;
//    }
//
//    public double getMinimumGoodAcceptanceLevel() {
//        return 0.10;
//    }
//
//    public double getMaximumGoodAcceptanceLevel() {
//        return 0.20;
//    }
//
//    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
//
//        public String getParserName() {
//            return UP_DOWN_OPERATOR;
//        }
//
//        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//
//            final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
//
//            final double weight = xo.getDoubleAttribute(WEIGHT);
//            CoercionMode mode = CoercionMode.parseMode(xo);
//
//            Parameter param1 = (Parameter) xo.getElementFirstChild(UP);
//            Parameter param2 = (Parameter) xo.getElementFirstChild(DOWN);
//
//            return new UpDownOperator(param1, param2, scaleFactor, weight, mode);
//        }
//
//        public String getParserDescription() {
//            return "This element represents an operator that scales two parameters in different directions. " +
//                    "Each operation involves selecting a scale uniformly at random between scaleFactor and 1/scaleFactor. " +
//                    "The up parameter is multipled by this scale and the down parameter is divided by this scale.";
//        }
//
//        public Class getReturnType() {
//            return UpDownOperator.class;
//        }
//
//        public XMLSyntaxRule[] getSyntaxRules() {
//            return rules;
//        }
//
//        private final XMLSyntaxRule[] rules = {
//                AttributeRule.newDoubleRule(SCALE_FACTOR),
//                AttributeRule.newDoubleRule(WEIGHT),
//                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
//                new ElementRule(UP,
//                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
//                new ElementRule(DOWN,
//                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
//        };
//    };
//
//    //PRIVATE STUFF
//
//    private Parameter upParameter = null;
//    private Parameter downParameter = null;
//    private final ContinuousVariablePrior prior = new ContinuousVariablePrior();
//    private double scaleFactor = 0.5;
//}
