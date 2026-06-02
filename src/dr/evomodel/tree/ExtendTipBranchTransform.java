/*
 * ExtendTipBranchTransform.java
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

package dr.evomodel.tree;

import dr.inference.model.Parameter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Variable;

public class ExtendTipBranchTransform extends TreeTransform {
    public ExtendTipBranchTransform(Parameter diffusionPrecision) {
        super("extendTipBranchTransform");
        diffusionPrecision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        this.diffusionPrecision = diffusionPrecision;
        addVariable(diffusionPrecision);
    }

    @Override
    public double transform(Tree tree, NodeRef node, double originalHeight) {
        double rootHeight = tree.getNodeHeight(tree.getRoot()) ;
        if (tree.isExternal(node)) {
            return originalHeight;
        }
        return originalHeight + getExtensionForNode(tree, node) * rootHeight;
    }

    @Override
    protected double getScaleForNode(Tree tree, NodeRef node) {
        return diffusionPrecision.getParameterValue(0);
    }

    private double getExtensionForNode(Tree tree, NodeRef node) {
        double samplingVariance = 1.0 - 1.0 / getScaleForNode(tree, node);
        return samplingVariance * getScaleForNode(tree, node);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(diffusionPrecision);
    }

    @Override
    public String getInfo() {
        return null;
    }

    private final Parameter diffusionPrecision;
}
