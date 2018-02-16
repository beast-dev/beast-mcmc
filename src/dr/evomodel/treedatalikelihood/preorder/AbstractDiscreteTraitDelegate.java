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

/**
 * AbstractDiscreteTraitDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class AbstractDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate {

    private final BeagleDataLikelihoodDelegate likelihoodDelegate;

    AbstractDiscreteTraitDelegate(String name,
                                  Tree tree,
                                  BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);
        this.likelihoodDelegate = likelihoodDelegate;
    }

    protected boolean isLoggable() {
        return true;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        // TODO
    }

    @Override
    public void modelRestored(Model model) {

    }

    public static String getName(String name) {
        return "derivative." + name;
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

    @Override
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

    protected double[] getTraitForNode(NodeRef node) {

        assert (node != null);

        // TODO
        return null;
    }

    @Override
    public void setupStatistics() {
        // TODO
    }

    @Override
    protected void simulateRoot(int rootNumber) {
        // TODO
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        // TODO
    }


    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        return likelihoodDelegate.vectorizeNodeOperations(nodeOperations, operations);
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
}

