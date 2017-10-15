/*
 * ScaleOperator.java
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

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.math.MathUtils;

import java.util.logging.Logger;

/**
 * A generic scale operator for use with a multi-dimensional parameters.
 * Either scale all dimentions at once or scale one dimention at a time.
 * An optional bit vector and a threshold is used to vary the rate of the individual dimentions according
 * to their on/off status. For example a threshold of 1 means pick only "on" dimentions.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ScaleOperator.java,v 1.20 2005/06/14 10:40:34 rambaut Exp $
 */
public class ScaleOperator extends AbstractCoercableOperator {

    private Parameter indicator;
    private double indicatorOnProb;

    public ScaleOperator(Variable variable, double scale) {

        this(variable, scale, CoercionMode.COERCION_ON, 1.0);
    }

    public ScaleOperator(Variable<Double> variable, double scale, CoercionMode mode, double weight) {

        this(variable, false, 0, scale, mode, null, 1.0, false);
        setWeight(weight);
    }

    public ScaleOperator(Variable<Double> variable, boolean scaleAll, int degreesOfFreedom, double scale,
                         CoercionMode mode, Parameter indicator, double indicatorOnProb, boolean scaleAllInd) {

        super(mode);

        this.variable = variable;
        this.indicator = indicator;
        this.indicatorOnProb = indicatorOnProb;
        this.scaleAll = scaleAll;
        this.scaleAllIndependently = scaleAllInd;
        this.scaleFactor = scale;
        this.degreesOfFreedom = degreesOfFreedom;
    }


    /**
     * @return the parameter this operator acts on.
     */
    public Variable getVariable() {
        return variable;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        double logq;

        final Bounds<Double> bounds = variable.getBounds();
        final int dim = variable.getSize();

        if (scaleAllIndependently) {
            // update all dimensions independently.
            logq = 0;
            for (int i = 0; i < dim; i++) {

                final double scaleOne = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
                final double offset = bounds.getLowerLimit(i);

                // scale offset by the lower bound
                final double value = ((variable.getValue(i) - offset) * scaleOne) + offset;

                logq -= Math.log(scaleOne);

                if (value > bounds.getUpperLimit(i)) {
                    // if bounded then perhaps this could be reflected.
                    throw new RuntimeException("proposed value greater than upper bound");
                }

                variable.setValue(i, value);

            }
        } else if (scaleAll) {
            // update all dimensions
            // hasting ratio is dim-2 times of 1dim case. would be nice to have a reference here
            // for the proof. It is supposed to be somewhere in an Alexei/Nicholes article.
            if (degreesOfFreedom > 0)
                // For parameters with non-uniform prior on only one dimension
                logq = -degreesOfFreedom * Math.log(scale);
            else
                logq = (dim - 2) * Math.log(scale);

            // Must first set all parameters first and check for boundaries later for the operator to work
            // correctly with dependent parameters such as tree node heights.
            for (int i = 0; i < dim; i++) {
                // For scale all we scale by the same factor (i.e., not relative to their individual
                // origins).
                variable.setValue(i, variable.getValue(i) * scale);
            }

            for (int i = 0; i < dim; i++) {
                if (variable.getValue(i) > variable.getBounds().getUpperLimit(i)) {
                    throw new RuntimeException("proposed value greater than upper bound");
                }
            }
        } else {
            logq = -Math.log(scale);

            // which bit to scale
            int index;
            if (indicator != null) {
                final int idim = indicator.getDimension();
                final boolean impliedOne = idim == (dim - 1);
                // available bit locations
                int[] loc = new int[idim + 1];
                int nLoc = 0;
                // choose active or non active ones?
                final boolean takeOne = indicatorOnProb >= 1.0 || MathUtils.nextDouble() < indicatorOnProb;

                if (impliedOne && takeOne) {
                    loc[nLoc] = 0;
                    ++nLoc;
                }
                for (int i = 0; i < idim; i++) {
                    final double value = indicator.getStatisticValue(i);
                    if (takeOne == (value > 0)) {
                        loc[nLoc] = i + (impliedOne ? 1 : 0);
                        ++nLoc;
                    }
                }

                if (nLoc > 0) {
                    final int rand = MathUtils.nextInt(nLoc);
                    index = loc[rand];
                } else {
                    // this used to throw an exception
                    return Double.NEGATIVE_INFINITY;
                }
            } else {
                // any is good
                index = MathUtils.nextInt(dim);
            }

            final double oldValue = variable.getValue(index);
            final double offset = bounds.getLowerLimit(index);

            if (oldValue == 0) {
                Logger.getLogger("dr.inference").severe("The " + ScaleOperatorParser.SCALE_OPERATOR +
                        " for " +
                        variable.getVariableName()
                        + " has failed since the parameter has a value of 0.0." +
                        "\nTo fix this problem, initalize the value of " +
                        variable.getVariableName() + " to be a positive real number"
                );
            }
            final double newValue = ((oldValue - offset) * scale) + offset;

            if (newValue > bounds.getUpperLimit(index)) {
                // if bounded then perhaps this could be reflected.
                throw new RuntimeException("proposed value greater than upper bound: " + newValue + " (" + variable.getId() + ")");
            }

            variable.setValue(index, newValue);

            // provides a hook for subclasses
            cleanupOperation(newValue, oldValue);
        }

        return logq;
    }

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     *
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "scale(" + variable.getVariableName() + ")";
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public String toString() {
        return ScaleOperatorParser.SCALE_OPERATOR + "(" + variable.getVariableName() + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private Variable<Double> variable = null;
    private boolean scaleAll = false;
    private boolean scaleAllIndependently = false;
    private int degreesOfFreedom = 0;
    private double scaleFactor = 0.5;
}
