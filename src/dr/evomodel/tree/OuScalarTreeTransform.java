/*
 * SingleScalarTreeTransform.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class OuScalarTreeTransform extends TreeTransform {

    public OuScalarTreeTransform(Parameter alpha) {
        super("OuScalarTreeTransform");
        alpha.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.alpha = alpha;
        addVariable(alpha);
    }

    public double transform(Tree tree, NodeRef node, double originalHeight) {
        if (getScaleForNode(tree, node)==0){
            return originalHeight;
        } else {
            if (tree.isExternal(node)) {
                return originalHeight;
            }
            double transformedHeight = (1 - Math.exp(-2*getScaleForNode(tree, node)*originalHeight))/(2*getScaleForNode(tree, node));
            final double rootHeight = tree.getNodeHeight(tree.getRoot());
            double transformedRootHeight = (1 - Math.exp(-2*getScaleForNode(tree, node)*rootHeight))/(2*getScaleForNode(tree, node));
            return transformedHeight*(rootHeight/transformedRootHeight);
        }
    }

    protected double getScaleForNode(Tree tree, NodeRef node) {
        return alpha.getParameterValue(0);
    }

    public String getInfo() {
        return "OU transform by " + alpha.getId();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(alpha);
    }

    private final Parameter alpha;
}
