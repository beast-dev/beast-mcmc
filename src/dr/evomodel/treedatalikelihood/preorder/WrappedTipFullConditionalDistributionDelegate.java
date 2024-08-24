/*
 * WrappedTipFullConditionalDistributionDelegate.java
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class WrappedTipFullConditionalDistributionDelegate extends
        AbstractFullConditionalDistributionDelegate {

    private static final String NAME_PREFIX = "wFcd";

    public WrappedTipFullConditionalDistributionDelegate(String name, Tree tree,
                                                         MultivariateDiffusionModel diffusionModel,
                                                         ContinuousTraitPartialsProvider dataModel,
                                                         ConjugateRootTraitPrior rootPrior,
                                                         ContinuousRateTransformation rateTransformation,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
    }

    public static String getName(String name) {
        return NAME_PREFIX + "." + name;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {

        TreeTrait<List<WrappedNormalSufficientStatistics>> baseTrait =
                new TreeTrait<List<WrappedNormalSufficientStatistics>>() {

            public String getTraitName() {
                return getName(name);
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() { return List.class; }

            public List<WrappedNormalSufficientStatistics> getTrait(Tree t, NodeRef node) {
                assert (tree == t);

                return getTraitForNode(node);
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formatted(getTrait(tree, node));
            }

            public boolean getLoggable() {
                return isLoggable();
            }
        };

        treeTraitHelper.addTrait(baseTrait);
    }

    private static String formatted(List<WrappedNormalSufficientStatistics> statistics) {

        StringBuilder sb = new StringBuilder();
        for (WrappedNormalSufficientStatistics stat : statistics) {
            sb.append(stat.toString()).append(";\t");
        }

        return sb.toString();
    }

    protected List<WrappedNormalSufficientStatistics> getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert dimPartial > 0;
        assert dimTrait > 0;

        simulationProcess.cacheSimulatedTraits(node);

        final int numberOfNodes = (node == null) ? tree.getExternalNodeCount() : 1;
        final int numberOfBuffers = (node == null) ? cdi.getBufferCount() : 1;
        double[] partial = new double[dimPartial * numTraits * numberOfBuffers];

        final int index = (node == null) ? -1 : likelihoodDelegate.getActiveNodeIndex(node.getNumber());
        cdi.getPreOrderPartial(index, partial);

        List<WrappedNormalSufficientStatistics> statistics = new ArrayList<WrappedNormalSufficientStatistics>();

        for (int n = 0; n < numberOfNodes; ++n) {

            int mapped = (node == null) ?
                    likelihoodDelegate.getActiveNodeIndex(tree.getExternalNode(n).getNumber()) : n;

            assert (numTraits == 1); // TODO Generalize
            
            statistics.add(new WrappedNormalSufficientStatistics(
                    partial, mapped, dimTrait, Pd, likelihoodDelegate.getPrecisionType())
            );
        }

        return statistics;
    }
}
