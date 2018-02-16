/*
 * DataSimulationDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.*;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.inference.model.Model;

import java.util.List;

import static dr.evolution.tree.TreeTrait.DA.factory;

/**
 * AbstractDiscreteTraitDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class AbstractDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate
        implements TreeTrait.TraitInfo<double[]> {

    private final BeagleDataLikelihoodDelegate likelihoodDelegate;

    AbstractDiscreteTraitDelegate(String name,
                                  Tree tree,
                                  BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);
        this.likelihoodDelegate = likelihoodDelegate;
    }

    @Override
    public void simulate(final int[] operations, final int operationCount,
                         final int rootNodeNumber) {

        super.simulate(operations, operationCount, rootNodeNumber); // TODO Should override this to compute pre-order partials

        // TODO
    }

    @Override
    public void setupStatistics() {
        throw new RuntimeException("Not used (?) with BEAGLE");
    }

    @Override
    protected void simulateRoot(int rootNumber) {
        throw new RuntimeException("Not used with BEAGLE");
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        throw new RuntimeException("Not used with BEAGLE");
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(factory(this));
    }


    @Override
    public String getTraitName() {
        return "derivative." + name;
    }

    @Override
    public TreeTrait.Intent getTraitIntent() {
        return TreeTrait.Intent.NODE;
    }

    @Override
    public Class getTraitClass() {
        return double[].class;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {

        assert (tree == this.tree);
        assert (node == null); // Implies: get trait for all nodes at same time

        // TODO See TipGradientViaFullConditionalDelegate.getTrait() as an example of using post- and pre-order partials together

        return null;
    }

    @Override
    public boolean isTraitLoggable() {
        return false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        // TODO
    }

    @Override
    public void modelRestored(Model model) {

    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        return likelihoodDelegate.vectorizeNodeOperations(nodeOperations, operations);
    }
}

