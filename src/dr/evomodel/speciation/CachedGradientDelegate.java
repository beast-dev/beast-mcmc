/*
 * CachedGradientDelegate.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.speciation;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

class CachedGradientDelegate extends AbstractModel implements TreeTrait<double[]> {

    private final SpeciationModelGradientProvider provider;
    private final BigFastTreeIntervals treeIntervals;
    private final SpeciationModel speciationModel;

    private double[] gradient;
    private double[] storedGradient;
    private boolean gradientKnown;
    private boolean storedGradientKnown;

    CachedGradientDelegate(EfficientSpeciationLikelihood likelihood) {
        super("cachedGradientDelegate");

        this.provider = likelihood.getGradientProvider();
        this.treeIntervals = likelihood.getTreeIntervals();
        this.speciationModel = likelihood.getSpeciationModel();

        addModel(this.treeIntervals);
        addModel(this.speciationModel);

        gradientKnown = false;
    }

    private double[] getGradientLogDensityImpl() {
        double[] gradient;
        if (speciationModel instanceof BirthDeathEpisodicSeriallySampledModel) {
            gradient = new double[5*((BirthDeathEpisodicSeriallySampledModel) speciationModel).numIntervals];
        } else {
            gradient = new double[5];
        }

        provider.precomputeGradientConstants(); // TODO hopefully get rid of this
        provider.updateModelValues(0);

        double[] modelBreakPoints = provider.getBreakPoints();
        assert modelBreakPoints[modelBreakPoints.length - 1] == Double.POSITIVE_INFINITY;

        int currentModelSegment = 0;

        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {

            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            final int nLineages = treeIntervals.getLineageCount(i);

            while (intervalEnd > modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?

                final double segmentIntervalEnd = modelBreakPoints[currentModelSegment];
                provider.processGradientModelSegmentBreakPoint(gradient, currentModelSegment, intervalStart, segmentIntervalEnd);
                intervalStart = segmentIntervalEnd;
                ++currentModelSegment;
                provider.updateModelValues(currentModelSegment);
            }

            // TODO Need to check for intervalStart == intervalEnd?
            // TODO Need to check for intervalStart == intervalEnd == 0.0?

            if (intervalEnd > intervalStart) {
//                System.err.println("interval: " + intervalStart + " -- " + intervalEnd);
                provider.processGradientInterval(gradient, currentModelSegment, intervalStart, intervalEnd, nLineages);
            }

            // Interval ends with a coalescent or sampling event at time intervalEnd
            if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {

                provider.processGradientSampling(gradient, currentModelSegment, intervalEnd);

            } else if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {

                provider.processGradientCoalescence(gradient, currentModelSegment, intervalEnd);

            } else {
                throw new RuntimeException("Birth-death tree includes non birth/death/sampling event.");
            }
        }

        // We've missed the first sample and need to add it back
        // TODO May we missed multiple samples @ t == 0.0?
        provider.processGradientSampling(gradient, 0, treeIntervals.getStartTime()); // TODO for-loop for models with multiple segments?

        // origin branch is a fake branch that doesn't exist in the tree, now compute its contribution
        provider.processGradientOrigin(gradient, currentModelSegment, treeIntervals.getTotalDuration());

        provider.logConditioningProbability(gradient);

        return gradient;
    }

    @Override
    public String getTraitName() {
        return EfficientSpeciationLikelihoodGradient.GRADIENT_KEY;
    }

    @Override
    public Intent getIntent() {
        return Intent.WHOLE_TREE;
    }

    @Override
    public Class getTraitClass() {
        return double[].class;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {
        if (!gradientKnown) {
            gradient = getGradientLogDensityImpl();
            gradientKnown = true;
        }
        return gradient;
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {
        return null;
    }

    @Override
    public boolean getLoggable() {
        return false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeIntervals || model == speciationModel) {
            gradientKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model: " + model.getId());
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new IllegalArgumentException("Unknown variable: " + variable.getId());
    }

    @Override
    protected void storeState() {
        if (gradient != null) {
            if (storedGradient == null) {
                storedGradient = new double[gradient.length];
            }
            System.arraycopy(gradient, 0, storedGradient, 0, gradient.length);
        }
        storedGradientKnown = gradientKnown;
    }

    @Override
    protected void restoreState() {
        double[] swap = gradient;
        gradient = storedGradient;
        storedGradient = swap;

        gradientKnown = storedGradientKnown;
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }
}
