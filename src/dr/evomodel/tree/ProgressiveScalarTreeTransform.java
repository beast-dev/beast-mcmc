/*
 * ProgressiveScalarTreeTransform.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A. Suchard
 */
public class ProgressiveScalarTreeTransform extends TreeTransform {

    public ProgressiveScalarTreeTransform(Parameter scale) {
        super("progressiveScalarTreeTransform");
        scale.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        this.scale = scale;
        addVariable(scale);
    }

    public double transform(TransformedTreeModel tree, NodeRef node, double originalHeight) {
        if (tree.isExternal(node)) {
            return originalHeight;
        }
        if (tree.isRoot(node)) {
            return tree.getRootHeightParameter().getParameterValue(0);
        }

        final double parentHeight = tree.getNodeHeight(tree.getParent(node));
        return parentHeight - scale.getParameterValue(0) * (parentHeight - originalHeight);
    }

    public String getInfo() {
        return "Linear, progressive transform by " + scale.getId();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(scale);
    }

    private final Parameter scale;

    // TODO Move to JUnitTest

    public static void main(String[] args) {
        try {
            NewickImporter importer = new NewickImporter("(0:3.0,(1:2.0,(2:1.0,3:1):1.0):1.0);");
            Tree tree = importer.importTree(null);

            Parameter scale = new Parameter.Default(0.5);

//            TreeTransform xform = new SingleScalarTreeTransform(scale);
            TreeTransform xform = new ProgressiveScalarTreeTransform(scale);

            TransformedTreeModel model = new TransformedTreeModel("tree", tree, xform);

            System.err.println(model.toString());

        } catch (Exception e) {

        }
    }
}
