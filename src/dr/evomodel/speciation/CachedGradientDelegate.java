/*
 * CachedGradientDelegate.java
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
import dr.util.Timer;
import dr.xml.Reportable;

class CachedGradientDelegate extends AbstractModel implements TreeTrait<double[]> {

    private final SpeciationModelGradientProvider provider;
    private final BigFastTreeIntervals treeIntervals;
    private final SpeciationModel speciationModel;

    public static final boolean MEASURE_RUN_TIME = false;
    public double gradientTime;

    public int gradientCounts;
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
        gradientTime = 0;
        gradientCounts = 0;
        gradientKnown = false;
    }

    private double[] getGradientLogDensityImpl() {
        Timer timer;
        if (MEASURE_RUN_TIME) {
            timer = new Timer();
            timer.start();
        }
        double[] gradient = new double[provider.getGradientLength()];

        provider.precomputeGradientConstants(); // TODO hopefully get rid of this
        provider.updateGradientModelValues(0);

        double[] modelBreakPoints = provider.getBreakPoints();
        assert modelBreakPoints[modelBreakPoints.length - 1] == Double.POSITIVE_INFINITY;

        int currentModelSegment = 0;

        while (treeIntervals.getStartTime() >= modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?
            ++currentModelSegment;
            speciationModel.updateLikelihoodModelValues(currentModelSegment);
        }

        provider.processGradientSampling(gradient, currentModelSegment, treeIntervals.getStartTime()); // TODO Fix for getStartTime() != 0.0

        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {

            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            final int nLineages = treeIntervals.getLineageCount(i);

            while (intervalEnd >= modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?

                final double segmentIntervalEnd = modelBreakPoints[currentModelSegment];
                provider.processGradientModelSegmentBreakPoint(gradient, currentModelSegment, intervalStart, segmentIntervalEnd, nLineages);
                intervalStart = segmentIntervalEnd;
                ++currentModelSegment;
                provider.updateGradientModelValues(currentModelSegment);
            }

            if (intervalEnd > intervalStart) {
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

        // origin branch is a fake branch that doesn't exist in the tree, now compute its contribution
        provider.processGradientOrigin(gradient, currentModelSegment, treeIntervals.getTotalDuration());

        provider.logConditioningProbability(currentModelSegment,gradient);

        if (MEASURE_RUN_TIME) {
            timer.stop();
            double timeInSeconds = timer.toNanoSeconds();
            gradientTime += timeInSeconds;
            gradientCounts += 1;
        }

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

    public double getGradientTime() {
        return gradientTime;
    }
}
