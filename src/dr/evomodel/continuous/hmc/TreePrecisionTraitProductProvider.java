/*
 * FullyConjugateTreeTipsPotentialDerivative.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public abstract class TreePrecisionTraitProductProvider extends AbstractModel
        implements PrecisionMatrixVectorProductProvider, Reportable {

    final Tree tree;
    final Parameter dataParameter;
    protected final int dimTrait;
    protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    TreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                             ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        super("treePrecisionTraitProductProvider");

        this.tree = treeDataLikelihood.getTree();
        this.dataParameter = likelihoodDelegate.getDataModel().getParameter();
        this.likelihoodDelegate = likelihoodDelegate;
        this.dimTrait = likelihoodDelegate.getTraitDim();

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }

    public Tree getTree() { return tree; }

    public ContinuousDataLikelihoodDelegate getLikelihoodDelegate() { return likelihoodDelegate; }

    public ContinuousTraitPartialsProvider getDataModel() { return likelihoodDelegate.getDataModel(); }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    @Override
    protected void storeState() {
        // Do nothing
    }

    @Override
    protected void restoreState() {
        // Do nothing
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    double[] expensiveProduct(Parameter vector, double[][] treeTraitPrecision) {

        int dim = treeTraitPrecision.length;
        double[] result = new double[dim];
        for (int row = 0; row < dim; ++row) {
            double sum = 0.0;
            for (int col = 0; col < treeTraitPrecision[row].length; ++col) {
                sum += treeTraitPrecision[row][col] * vector.getParameterValue(col);
            }
            result[row] = sum;
        }
        return result;
    }

    void debug(double[] result, Parameter vector) {
        double[] result2 = expensiveProduct(vector, likelihoodDelegate.getTreeTraitPrecision());
        System.err.println("via FCD: " + new WrappedVector.Raw(result));
        System.err.println("direct : " + new WrappedVector.Raw(result2));
        System.err.println();
    }

    @Override
    public String getReport() {
        double[] result = getProduct(dataParameter);
        return (new WrappedVector.Raw(result, 0, result.length)).toString();
    }
}
