/*
 * ColouredOperator.java
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

package dr.evomodel.operators;

import dr.evomodel.coalescent.structure.ColourSamplerModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * An operator that wraps another operator and performs a re-colouring of the tree
 * immediately afterwards.
 *
 * @author Alexei Drummond
 * @version $Id: ColouredOperator.java,v 1.10 2006/07/28 11:41:23 rambaut Exp $
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class ColouredOperator implements CoercableMCMCOperator {

    public static final double ACCEPTANCE_FACTOR = 0.5;

    public static final String COLOURED_OPERATOR = "colouredOperator";

    private final ColourSamplerModel colouringModel;

    private final MCMCOperator innerOperator;

    public ColouredOperator(ColourSamplerModel colouringModel, MCMCOperator operator) {

        this.colouringModel = colouringModel;
        this.innerOperator = operator;
    }

    public final double operate() {

        double logP = colouringModel.getTreeColouring().getLogProbabilityDensity();

        double logO = innerOperator.operate();

        colouringModel.resample();
        double logQ = colouringModel.getTreeColouring().getLogProbabilityDensity();

        return logO + logP - logQ;
    }

    public double getCoercableParameter() {
        if (innerOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) innerOperator).getCoercableParameter();
        }
        throw new IllegalArgumentException();
    }

    public void setCoercableParameter(double value) {
        if (innerOperator instanceof CoercableMCMCOperator) {
            ((CoercableMCMCOperator) innerOperator).setCoercableParameter(value);
            return;
        }
        throw new IllegalArgumentException();
    }

    public double getRawParameter() {

        if (innerOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) innerOperator).getRawParameter();
        }
        throw new IllegalArgumentException();
    }

    public CoercionMode getMode() {
        if (innerOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) innerOperator).getMode();
        }
        return CoercionMode.COERCION_OFF;
    }

    public String getOperatorName() {
        return "Coloured(" + innerOperator.getOperatorName() + ")";
    }

    public Element createOperatorElement(Document d) {
        throw new RuntimeException("not implemented");
    }

//    public double getTargetAcceptanceProbability() { return 0.117; }
//	public double getMinimumAcceptanceLevel() { return 0.01; }
//	public double getMinimumGoodAcceptanceLevel() { return 0.05; }
//    public double getMaximumAcceptanceLevel() { return 0.40; }
//    public double getMaximumGoodAcceptanceLevel() { return 0.30; }

    public double getTargetAcceptanceProbability() {
        return innerOperator.getTargetAcceptanceProbability() * ACCEPTANCE_FACTOR;
    }

    public double getMinimumAcceptanceLevel() {
        return innerOperator.getMinimumAcceptanceLevel() * ACCEPTANCE_FACTOR;
    }

    public double getMaximumAcceptanceLevel() {
        return innerOperator.getMaximumAcceptanceLevel() * ACCEPTANCE_FACTOR;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return innerOperator.getMinimumGoodAcceptanceLevel() * ACCEPTANCE_FACTOR;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return innerOperator.getMaximumGoodAcceptanceLevel() * ACCEPTANCE_FACTOR;
    }

    // All of this is copied and modified from SimpleMCMCOperator
    /**
     * @return the weight of this operator.
     */
    public final double getWeight() {
        return innerOperator.getWeight();
    }

    /**
     * Sets the weight of this operator.
     */
    public final void setWeight(double w) {
        innerOperator.setWeight(w);
    }

    public final void accept(double deviation) {
        innerOperator.accept(deviation);
    }

    public final void reject() {
        innerOperator.reject();
    }

    public final void reset() {
        innerOperator.reset();
    }

    public final long getCount() {
        return innerOperator.getCount();
    }

    public final long getAcceptCount() {
        return innerOperator.getAcceptCount();
    }

    public final void setAcceptCount(long accepted) {
        innerOperator.setAcceptCount(accepted);
    }

    public final long getRejectCount() {
        return innerOperator.getRejectCount();
    }

    public final void setRejectCount(long rejected) {
        innerOperator.setRejectCount(rejected);
    }

    public final double getMeanDeviation() {
        return innerOperator.getMeanDeviation();
    }

    public final double getSumDeviation() {
        return innerOperator.getSumDeviation();
    }

    public double getSpan(boolean reset) {
        return 0;
    }

    public final void setSumDeviation(double sumDeviation) {
        innerOperator.setSumDeviation(sumDeviation);
    }

    public String getPerformanceSuggestion() {
        return innerOperator.getPerformanceSuggestion();
    }

    public double getMeanEvaluationTime() {
        return innerOperator.getMeanEvaluationTime();
    }

    public long getTotalEvaluationTime() {
        return innerOperator.getTotalEvaluationTime();
    }

    public void addEvaluationTime(long time) {
        innerOperator.addEvaluationTime(time);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COLOURED_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) {

            MCMCOperator operator = (MCMCOperator) xo.getChild(MCMCOperator.class);
            ColourSamplerModel colourSamplerModel = (ColourSamplerModel) xo.getChild(ColourSamplerModel.class);

            return new ColouredOperator(colourSamplerModel, operator);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an arbitrary operator immediately followed by a re-colouring of the tree";
        }

        public Class getReturnType() {
            return ColouredOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MCMCOperator.class),
                new ElementRule(ColourSamplerModel.class)
        };

    };
}
