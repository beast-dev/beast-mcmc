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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;

import java.util.Set;

/**
 * @author Andy Magee
 * @author Yucai Shao
 * @author Marc Suchard
 */
public class EfficientSpeciationLikelihood extends SpeciationLikelihood implements TreeTraitProvider {

    private final BigFastTreeIntervals treeIntervals;
    private final TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper();

    private boolean intervalsKnown;

    private final double TOLERANCE = 1e-5;

    public EfficientSpeciationLikelihood(Tree tree, SpeciationModel speciationModel, Set<Taxon> exclude, String id) {
        super(tree, speciationModel, exclude, id);

        if (!(tree instanceof TreeModel)) {
            throw new IllegalArgumentException("Must currently provide a TreeModel");
        }

//        fixTimes();

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

        logL += speciationModel.logConditioningProbability();

        return logL;
    }

    private void fixTimes() {

        FlexibleTree binaryTree = new FlexibleTree(tree, true);

        double[] intervalTimes = speciationModel.getBreakPoints();
        for (int i = 0; i < binaryTree.getExternalNodeCount(); i++) {
            // TODO we can be lazy since we only do this once but a linear search is still sad
            NodeRef node = binaryTree.getNode(i);
            double thisTipTime = binaryTree.getNodeHeight(node);
//            System.err.println("Working on tip " + i + " at time " + thisTipTime);
            // TODO
            if (thisTipTime < TOLERANCE) {
//                System.err.println("Adusting time " + thisTipTime + " to 0.0");
                binaryTree.setNodeHeight(node,0.0);
            } else {
                for (int j = 0; j < intervalTimes.length; j++) {
                    if (Math.abs(thisTipTime - intervalTimes[j]) < TOLERANCE) {
//                        System.err.println("Adusting time " + thisTipTime + " to " + intervalTimes[j]);
                        binaryTree.setNodeHeight(node,intervalTimes[j]);
                        break;
                    }
                }
            }
        }
        tree = binaryTree;
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
}