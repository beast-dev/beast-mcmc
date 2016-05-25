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
 * @author Marc A. Suchard
 */
public class SingleScalarTreeTransform extends TreeTransform {

    public SingleScalarTreeTransform(Parameter scale) {
        super("singleScalarTreeTransform");
        scale.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        this.scale = scale;
        addVariable(scale);
    }

    public double transform(Tree tree, NodeRef node, double originalHeight) {
        if (tree.isExternal(node)) {
            return originalHeight;
        }
        final double rootHeight = tree.getNodeHeight(tree.getRoot());
        return rootHeight - getScaleForNode(tree, node) * (rootHeight - originalHeight);
    }

    protected double getScaleForNode(Tree tree, NodeRef node) {
        return scale.getParameterValue(0);
    }

    public String getInfo() {
        return "Linear transform by " + scale.getId();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(scale);
    }

    private final Parameter scale;
}
