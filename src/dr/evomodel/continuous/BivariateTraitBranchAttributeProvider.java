/*
 * BivariateTraitBranchAttributeProvider.java
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

package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;

/**
 * @author Marc Suchard
 */
public abstract class BivariateTraitBranchAttributeProvider extends TreeTrait.DefaultBehavior implements TreeTrait<Double> {

    public static final String FORMAT = "%5.4f";

    public BivariateTraitBranchAttributeProvider(AbstractMultivariateTraitLikelihood traitLikelihood) {

        traitName = traitLikelihood.getTraitParameter().getId();
        label = traitName + extensionName();

        double[] rootTrait = traitLikelihood.getRootNodeTrait();
        if (rootTrait.length != 2)
            throw new RuntimeException("BivariateTraitBranchAttributeProvider only works for 2D traits");
    }

    protected abstract String extensionName();

    protected double branchFunction(double[] startValue, double[] endValue, double startTime, double endTime) {
        return convert(endValue[0]-startValue[0], endValue[1] - startValue[1], startTime - endTime);
    }

    protected abstract double convert(double latDifference, double longDifference, double timeDifference);

    public String getTraitName() {
        return label;
    }

    public Intent getIntent() {
        return Intent.BRANCH;
    }

    public Class getTraitClass() {
        return Double.class;
    }

    public int getDimension() {
        return 1;
    }

    public Double getTrait(Tree tree, NodeRef node) {
        if (tree != traitLikelihood.getTreeModel())
            throw new RuntimeException("Bad bug.");

        NodeRef parent = tree.getParent(node);
        double[] startTrait = traitLikelihood.getTraitForNode(tree, parent, traitName);
        double[] endTrait = traitLikelihood.getTraitForNode(tree, node, traitName);
        double startTime = tree.getNodeHeight(parent);
        double endTime = tree.getNodeHeight(node);

        return branchFunction(startTrait, endTrait, startTime, endTime);
    }

    public String getTraitString(Tree tree, NodeRef node) {
        NodeRef parent = tree.getParent(node);
        double[] startTrait = traitLikelihood.getTraitForNode(tree, parent, traitName);
        double[] endTrait = traitLikelihood.getTraitForNode(tree, node, traitName);
        double startTime = tree.getNodeHeight(parent);
        double endTime = tree.getNodeHeight(node);

        return String.format(BivariateTraitBranchAttributeProvider.FORMAT,
                branchFunction(startTrait, endTrait, startTime, endTime));
    }

    protected AbstractMultivariateTraitLikelihood traitLikelihood;
    protected String traitName;
    protected String label;


}
