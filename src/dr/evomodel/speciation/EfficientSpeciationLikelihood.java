/*
 * EfficientSpeciationLikelihood.java
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
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.bigfasttree.ModelCompressedBigFastTreeIntervals;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.math.MathUtils;
import dr.math.distributions.BetaDistribution;
import dr.util.Timer;

import java.util.Set;

/**
 * @author Andy Magee
 * @author Yucai Shao
 * @author Marc Suchard
 */
public class EfficientSpeciationLikelihood extends SpeciationLikelihood implements TreeTraitProvider {

    private final BigFastTreeIntervals treeIntervals;
    private final TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper();

    public static final boolean MEASURE_RUN_TIME = false;
    public double likelihoodTime;
    public int likelihoodCounts;

    private boolean intervalsKnown;

    private final double TOLERANCE = 1e-5;

    public EfficientSpeciationLikelihood(Tree tree, SpeciationModel speciationModel, Set<Taxon> exclude, String id) {
        super(tree, speciationModel, exclude, id);

        if (!(tree instanceof DefaultTreeModel)) {
            throw new IllegalArgumentException("Must currently provide a DefaultTreeModel");
        }

        likelihoodTime = 0;
        likelihoodCounts = 0;

        fixTimes();

        treeIntervals = new BigFastTreeIntervals((TreeModel)tree);

        addModel(treeIntervals);
    }

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);
        if (model == treeIntervals) {
            intervalsKnown = false;
        }
//        fireModelChanged(object, index);
    }

    final TreeModel getTreeModel() {
        return (TreeModel) tree;
    }

    final BigFastTreeIntervals getTreeIntervals() {
        return treeIntervals;
    }

    @Override
    double calculateLogLikelihood() {
        Timer timer;
        if (MEASURE_RUN_TIME) {
            timer = new Timer();
            timer.start();
        }
        speciationModel.updateLikelihoodModelValues(0);

        double[] modelBreakPoints = speciationModel.getBreakPoints();
        assert modelBreakPoints[modelBreakPoints.length - 1] == Double.POSITIVE_INFINITY;

        int currentModelSegment = 0;

        double logL = speciationModel.processSampling(0, treeIntervals.getStartTime()); // TODO Fix for getStartTime() != 0.0

        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {

            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            final int nLineages = treeIntervals.getLineageCount(i);

            while (intervalEnd >= modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?

                final double segmentIntervalEnd = modelBreakPoints[currentModelSegment];
                logL += speciationModel.processModelSegmentBreakPoint(currentModelSegment, intervalStart, segmentIntervalEnd, nLineages);
                intervalStart = segmentIntervalEnd;
                ++currentModelSegment;
                speciationModel.updateLikelihoodModelValues(currentModelSegment);
            }

            if (intervalEnd > intervalStart) {
                logL += speciationModel.processInterval(currentModelSegment, intervalStart, intervalEnd, nLineages);
            }

            // Interval ends with a coalescent or sampling event at time intervalEnd
            if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {

                logL += speciationModel.processSampling(currentModelSegment, intervalEnd);

            } else if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {

                logL += speciationModel.processCoalescence(currentModelSegment,intervalEnd);

            } else {
                throw new RuntimeException("Birth-death tree includes non birth/death/sampling event.");
            }
        }

        // origin branch is a fake branch that doesn't exist in the tree, now compute its contribution
        logL += speciationModel.processOrigin(currentModelSegment, treeIntervals.getTotalDuration());

        logL += speciationModel.logConditioningProbability(currentModelSegment);

        if (MEASURE_RUN_TIME) {
            timer.stop();
            double timeInSeconds = timer.toNanoSeconds();
            likelihoodTime += timeInSeconds;
            likelihoodCounts += 1;
        }

        return logL;
    }

    private void fixTimes() {

        // DefaultTreeModel cleanTree = new DefaultTreeModel(tree);
        TreeModel cleanTree = getTreeModel();

        double[] intervalTimes = speciationModel.getBreakPoints();
        for (int i = 0; i < cleanTree.getExternalNodeCount(); i++) {
            // TODO we can be lazy since we only do this once but a linear search is still sad
            NodeRef node = cleanTree.getExternalNode(i);
            double thisTipTime = cleanTree.getNodeHeight(node);
//            System.err.println("Working on tip " + i + " at time " + thisTipTime);
            // TODO
            if (thisTipTime < TOLERANCE) {
//                System.err.println("Adusting time " + thisTipTime + " to 0.0");
                cleanTree.setNodeHeight(node,0.0);
            } else {
                for (int j = 0; j < intervalTimes.length; j++) {
                    if (Math.abs(thisTipTime - intervalTimes[j]) < TOLERANCE) {
//                        System.err.println("Adusting time " + thisTipTime + " to " + intervalTimes[j]);
                        cleanTree.setNodeHeight(node,intervalTimes[j]);
                        break;
                    }
                }
            }
        }
        boolean adjustedBirths = false;
        for (int i = 0; i < cleanTree.getInternalNodeCount(); i++) {
            // TODO we can be lazy since we only do this once but a linear search is still sad
            NodeRef node = cleanTree.getInternalNode(i);
            double thisNodeTime = cleanTree.getNodeHeight(node);
            if (thisNodeTime == 0.0) {
                adjustedBirths = true;
                System.err.println("Some births were found at time 0.0 and moved (by no more than " + TOLERANCE + ") to avoid numerical issues.");
                cleanTree.setNodeHeight(node, TOLERANCE * MathUtils.nextBeta(2.0,1.0));
            } else {
                for (int j = 0; j < intervalTimes.length; j++) {
                    if (thisNodeTime == intervalTimes[j]) {
                        adjustedBirths = true;
                        double bound = thisNodeTime > TOLERANCE ? TOLERANCE : thisNodeTime;
                        double dt = bound * MathUtils.nextBeta(2.0,1.0);
                        if( MathUtils.nextBoolean() ) {
                            dt = -dt;
                        }
//                        System.err.println("Adusting time " + thisNodeTime + " to " + (thisNodeTime + dt));
                        cleanTree.setNodeHeight(node,thisNodeTime + dt);
                        break;
                    }
                }
            }
        }
        if (adjustedBirths) {
            System.err.println("Some births were exactly at event-sampling times and have been moved (by no more than " + TOLERANCE + ") to avoid numerical issues.");
        }
        // tree = cleanTree;
//        System.err.println("Adjusted tip times to match interval times.");
    }

    // Super-clean interface (just one intrusive function) and a better place, since `Likelihood`s have gradients (`Model`s do not).
    public SpeciationModelGradientProvider getGradientProvider() {
        if (gradientProvider == null) {
            gradientProvider = speciationModel.getProvider();
        }
        return gradientProvider;
    }

    private SpeciationModelGradientProvider gradientProvider = null;

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public void addTrait(TreeTrait trait) {
        treeTraits.addTrait(trait);
    }

    public String getReport() {
        String message = super.getReport();
        if (MEASURE_RUN_TIME) {
            message += "\n";
            // add likelihood calculation time
            message += "Likelihood calculation time is " + likelihoodTime / likelihoodCounts + " nanoseconds.\n";
        }
        return message;
    }
}