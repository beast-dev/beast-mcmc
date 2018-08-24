/*
 * ProgressiveScalarTreeTransform.java
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

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class ProgressiveScalarTreeTransform extends TreeTransform {

    public ProgressiveScalarTreeTransform(Parameter scale) {
        this(scale, null);
    }

    public ProgressiveScalarTreeTransform(MutableTreeModel tree, Parameter scale) {
        this(scale, new TrialTreeParameterModel(tree, scale, false, false, TreeTrait.Intent.BRANCH));
    }

    public ProgressiveScalarTreeTransform(Parameter scale, TrialTreeParameterModel treeParameterModel) {
        super("progressiveScalarTreeTransform");
        int dim = 1;
        this.treeParameterModel = treeParameterModel;
        this.scale = scale;

        if (treeParameterModel != null) {
            dim = treeParameterModel.getParameterSize();
            addModel(treeParameterModel);
        } else {
            addVariable(scale);
        }
        scale.addBounds(new Parameter.DefaultBounds(1.0, 0.0, dim));

        StringBuilder sb = new StringBuilder("Creating a branch-specific phenotypic mixture model.\n");
        sb.append("\tPlease cite:\n");
        sb.append(Citable.Utils.getCitationString(this));

        java.util.logging.Logger.getLogger("dr.evomodel.tree").info(sb.toString());
    }

    public double transform(Tree tree, NodeRef node, double originalHeight) {
        // Early exit
        if (tree.isExternal(node)) {
            return originalHeight;
        }
        if (tree.isRoot(node)) {
            return tree.getNodeHeight(node);
        }

        // Do recursive work
        final double parentHeight = tree.getNodeHeight(tree.getParent(node));
        return parentHeight - getScaleForNode(tree, node) * (parentHeight - originalHeight);
    }

    protected double getScaleForNode(Tree tree, NodeRef node) {
        if (treeParameterModel != null) {
            return treeParameterModel.getNodeValue(tree, node);
        } else {
            return scale.getParameterValue(0);
        }
    }

    public String getInfo() {
        return "Linear, progressive transform by " + scale.getId();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(treeParameterModel);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(scale);
    }

    private final Parameter scale;
    private final TrialTreeParameterModel treeParameterModel;

    // TODO Move to JUnitTest
    public static void main(String[] args) throws Exception {

        NewickImporter importer = new NewickImporter(
                "((((A:2,B:1):0.5,C:3):1.5,(D:1.5,E:1):3.5):1,F:5);"
//                "(0:3.0,(1:2.0,(2:1.0,3:1):1.0):1.0);"
        );
        Tree tree = importer.importTree(null);
        Parameter scale = new Parameter.Default(0.5);
        TreeTransform xform = new ProgressiveScalarTreeTransform(scale);
        TreeModel treeModel = new TreeModel("original", tree);
        TransformedTreeModel model = new TransformedTreeModel("tree", treeModel, xform);
        System.err.println(model.toString());


        TreeTransform xform2 = new SingleScalarTreeTransform(scale);
        TransformedTreeModel model2 = new TransformedTreeModel("tree2", treeModel, xform2);
        System.err.println(model2.toString());
    }
}
