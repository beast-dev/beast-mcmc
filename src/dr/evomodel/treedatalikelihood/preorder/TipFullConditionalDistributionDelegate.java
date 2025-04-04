/*
 * TipFullConditionalDistributionDelegate.java
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

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class TipFullConditionalDistributionDelegate extends AbstractFullConditionalDistributionDelegate {

    public TipFullConditionalDistributionDelegate(String name,
                                                  Tree tree,
                                                  MultivariateDiffusionModel diffusionModel,
                                                  ContinuousTraitPartialsProvider dataModel,
                                                  ConjugateRootTraitPrior rootPrior,
                                                  ContinuousRateTransformation rateTransformation,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
    }

    public int vectorizeNodeOperations(final List<NodeOperation> nodeOperations,
                                       final int[] operations) {
        return likelihoodDelegate.vectorizeNodeOperations(nodeOperations, operations);
    }

    public static String getName(String name) {
        return "fcd." + name;
    }

    public String getTraitName(String name) {
        return getName(name);
    }

    private String delegateGetTraitName() {
        return getTraitName(name);
    }

    private Class delegateGetTraitClass() {
        return double[].class;
    }

    protected void constructTraits(Helper treeTraitHelper) {

        TreeTrait.DA baseTrait = new TreeTrait.DA() {

            public String getTraitName() {
                return delegateGetTraitName();
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() {
                return delegateGetTraitClass();
            }

            public double[] getTrait(Tree t, NodeRef node) {
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

    private static String formatted(double[] values) {

        if (values.length == 1) {
            return Double.toString(values[0]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("{");

        for (int i = 0; i < values.length; ++i) {
            sb.append(Double.toString(values[i]));
            if (i < (values.length - 1)) {
                sb.append(",");
            }
        }

        sb.append("}");
        return sb.toString();
    }

    protected double[] getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert node != null;

        simulationProcess.cacheSimulatedTraits(node);

        double[] partial = new double[dimPartial * numTraits];
        cdi.getPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), partial);

        return partial;
    }
}
