/*
 * SpeciationLikelihood.java
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

package dr.evomodel.speciation;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;

import java.util.Set;

/**
 * @author Andy Magee
 * @author Yucai Shao
 * @author Marc Suchard
 */
public class EfficientSpeciationLikelihood extends SpeciationLikelihood {

    public EfficientSpeciationLikelihood(Tree tree, SpeciationModel speciationModel, Set<Taxon> exclude, String id) {
        super(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, tree, speciationModel, exclude);
        setId(id);
    }

    @Override
    double calculateLogLikelihood() {

        if (!(tree instanceof TreeModel)) {
            throw new IllegalArgumentException("Failed test");
        }

        speciationModel.precomputeConstants();

        // TODO Make cached class-object
        BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals((TreeModel)tree);

        double[] modelBreakPoints = speciationModel.getBreakPoints();
        assert modelBreakPoints[modelBreakPoints.length - 1] == Double.POSITIVE_INFINITY;

        int currentModelSegment = 0;

        double logL = 0.0;

        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {

            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            final int nLineages = treeIntervals.getLineageCount(i);

            while (intervalEnd > modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?

                final double segmentIntervalEnd = modelBreakPoints[currentModelSegment];
                logL += speciationModel.processModelSegmentBreakPoint(currentModelSegment, intervalStart, segmentIntervalEnd);
                intervalStart = segmentIntervalEnd;
                ++currentModelSegment;
            }

            // TODO Need to check for intervalStart == intervalEnd?
            // TODO Need to check for intervalStart == intervalEnd == 0.0?

            logL += speciationModel.processInterval(currentModelSegment, intervalStart, intervalEnd, nLineages);

            // Interval ends with a coalescent or sampling event at time intervalEnd
            if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {

                logL += speciationModel.processSampling(currentModelSegment, intervalEnd);

            } else if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {

                logL += speciationModel.processCoalescence(currentModelSegment,intervalEnd);

            } else {
                throw new RuntimeException("Birth-death tree includes non birth/death/sampling event.");
            }
        }

        // We've missed the first sample and need to add it back
        // TODO May we missed multiple samples @ t == 0.0?
        logL += speciationModel.processSampling(0, treeIntervals.getStartTime()); // TODO for-loop for models with multiple segments?

        // origin branch is a fake branch that doesn't exist in the tree, now compute its contribution
        logL += speciationModel.processOrigin(currentModelSegment, treeIntervals.getTotalDuration());

        logL += speciationModel.logConditioningProbability();

        return logL;
    }

    // Super-clean interface (just one intrusive function) and a better place, since `Likelihood`s have gradients (`Model`s do not).
    public SpeciationModelGradientProvider getGradientProvider() {
        if (gradientProvider == null) {
            gradientProvider = speciationModel.getProvider();
        }
        return gradientProvider;
    }

    private SpeciationModelGradientProvider gradientProvider = null;
}