/*
 * EmbeddedOrdinalLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.distribution;

import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.hmc.EmbeddedOrdinalParameter;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Likelihood for one or more embedded ordinal latent parameters sharing a
 * common set of interval cuts and ordinal log weights.
 *
 * Each continuous latent coordinate is mapped to an ordinal state by interval
 * cuts, and each state is assigned a log weight. The total log density is the
 * sum of the coordinatewise ordinal contributions. This gives a simple but
 * genuinely multi-dimensional discontinuous target with exact local updates.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class EmbeddedOrdinalLikelihood extends AbstractModelLikelihood implements DiscontinuousPotentialProvider {

    public static final String EMBEDDED_ORDINAL_LIKELIHOOD = "embeddedOrdinalLikelihood";

    private final Parameter latentParameter;
    private final Parameter logWeightParameter;
    private final Parameter cutParameter;
    private EmbeddedOrdinalParameter embedding;

    private double[] contributions;
    private double finiteSum;
    private int infiniteCount;

    private double[] storedContributions;
    private double storedFiniteSum;
    private int storedInfiniteCount;

    public EmbeddedOrdinalLikelihood(Parameter latentParameter,
                                     Parameter logWeightParameter,
                                     Parameter cutParameter) {
        super(EMBEDDED_ORDINAL_LIKELIHOOD);

        if (cutParameter.getDimension() != logWeightParameter.getDimension() + 1) {
            throw new IllegalArgumentException("The number of cuts must be one greater than the number of log weights");
        }

        this.latentParameter = latentParameter;
        this.logWeightParameter = logWeightParameter;
        this.cutParameter = cutParameter;

        addVariable(latentParameter);
        addVariable(logWeightParameter);
        addVariable(cutParameter);

        refreshEmbedding();
        recomputeContributions();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No child models.
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == cutParameter) {
            refreshEmbedding();
            recomputeContributions();
        } else if (variable == logWeightParameter) {
            recomputeContributions();
        } else if (variable == latentParameter) {
            if (type == Parameter.ChangeType.VALUE_CHANGED && index >= 0) {
                updateContribution(index);
            } else {
                recomputeContributions();
            }
        }
    }

    @Override
    protected void storeState() {
        storedContributions = contributions.clone();
        storedFiniteSum = finiteSum;
        storedInfiniteCount = infiniteCount;
    }

    @Override
    protected void restoreState() {
        final double[] temp = contributions;
        contributions = storedContributions;
        storedContributions = temp;
        finiteSum = storedFiniteSum;
        infiniteCount = storedInfiniteCount;
    }

    @Override
    protected void acceptState() {
        // contributions/finiteSum/infiniteCount already reflect the accepted state.
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        return getLogDensity();
    }

    @Override
    public void makeDirty() {
        recomputeContributions();
    }

    @Override
    public Parameter getParameter() {
        return latentParameter;
    }

    @Override
    public int getDimension() {
        return latentParameter.getDimension();
    }

    @Override
    public double getLogDensity() {
        return infiniteCount > 0 ? Double.NEGATIVE_INFINITY : finiteSum;
    }

    @Override
    public double getLogDensityAfterSingleCoordinateMove(int index, double proposedValue) {
        if (index < 0 || index >= latentParameter.getDimension()) {
            throw new IllegalArgumentException("Coordinate index out of bounds: " + index);
        }

        final double oldContribution = contributions[index];
        final boolean oldWasInfinite = isNegativeInfinity(oldContribution);
        final double newContribution = getLogDensityForCoordinate(proposedValue);

        final int otherInfiniteCount = infiniteCount - (oldWasInfinite ? 1 : 0);
        if (otherInfiniteCount > 0 || isNegativeInfinity(newContribution)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double finiteSumExcludingIndex = oldWasInfinite ? finiteSum : finiteSum - oldContribution;
        return finiteSumExcludingIndex + newContribution;
    }

    @Override
    public double getNextDiscontinuity(int index, double currentValue, double direction) {
        if (index < 0 || index >= latentParameter.getDimension()) {
            throw new IllegalArgumentException("Coordinate index out of bounds: " + index);
        }
        return embedding.getNextBoundary(currentValue, direction);
    }

    public int getOrdinalState() {
        if (latentParameter.getDimension() != 1) {
            throw new IllegalStateException("Use getOrdinalState(int) for multi-dimensional embedded ordinal likelihoods");
        }
        return getOrdinalState(0);
    }

    public int getOrdinalState(int index) {
        if (index < 0 || index >= latentParameter.getDimension()) {
            throw new IllegalArgumentException("Coordinate index out of bounds: " + index);
        }
        return embedding.getStateIndex(latentParameter.getParameterValue(index));
    }

    private double getLogDensityForCoordinate(double latentValue) {
        try {
            int state = embedding.getStateIndex(latentValue);
            return logWeightParameter.getParameterValue(state);
        } catch (IllegalArgumentException e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    private void refreshEmbedding() {
        embedding = new EmbeddedOrdinalParameter(cutParameter.getParameterValues());
    }

    private void recomputeContributions() {
        final int dimension = latentParameter.getDimension();
        contributions = new double[dimension];
        finiteSum = 0.0;
        infiniteCount = 0;
        for (int i = 0; i < dimension; i++) {
            setContribution(i, getLogDensityForCoordinate(latentParameter.getParameterValue(i)));
        }
    }

    private void updateContribution(int index) {
        final double oldContribution = contributions[index];
        if (isNegativeInfinity(oldContribution)) {
            infiniteCount--;
        } else {
            finiteSum -= oldContribution;
        }
        setContribution(index, getLogDensityForCoordinate(latentParameter.getParameterValue(index)));
    }

    private void setContribution(int index, double value) {
        contributions[index] = value;
        if (isNegativeInfinity(value)) {
            infiniteCount++;
        } else {
            finiteSum += value;
        }
    }

    private static boolean isNegativeInfinity(double value) {
        return Double.isInfinite(value) && value < 0.0;
    }
}
