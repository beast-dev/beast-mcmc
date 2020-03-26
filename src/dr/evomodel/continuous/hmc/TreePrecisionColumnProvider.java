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
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.ElementaryVectorDataModel;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import java.util.HashMap;
import java.util.Map;

import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public class TreePrecisionColumnProvider extends AbstractModel
        implements PrecisionColumnProvider, Reportable {

    private final TreePrecisionTraitProductProvider productProvider;

    final Tree tree;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final ContinuousTraitPartialsProvider tipData;
    private final Map<Integer, double[]> treeCache = new HashMap<>();

    private final int numTaxa;
    private final int dimTrait;

    public TreePrecisionColumnProvider(TreePrecisionTraitProductProvider productProvider) {

        super("treePrecisionColumnProvider");

        this.productProvider = productProvider;

        this.tree = productProvider.getTree();
        this.likelihoodDelegate = productProvider.likelihoodDelegate;
        this.tipData = productProvider.getDataModel();

        this.numTaxa = tree.getExternalNodeCount();
        this.dimTrait = likelihoodDelegate.getTraitDim();

        assert (likelihoodDelegate.getTraitCount() == 1);

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }

    @Override
    public double[] getColumn(int index) {

        double[] column = treeCache.get(index);
        if (column == null) {

            column = setDataModelAndGetColumn(index);
            treeCache.put(index, column);

            if (DEBUG_CACHE) {
                System.err.println("treePrecision key " + index + " not found " + new WrappedVector.Raw(column));
            }

        } else {
            if (DEBUG_CACHE) {
                System.err.println("treePrecision key " + index + " found    " + new WrappedVector.Raw(column));
            }
        }

        if (DEBUG) {
            debug(column, index);
        }

        return column;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            treeCache.clear();
        }
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

    private double[] setDataModelAndGetColumn(int index) {

        Parameter parameter = tipData.getParameter();
        WrappedVector savedParameter;

        if (RESET_DATA) {
            savedParameter = new WrappedVector.Raw(parameter.getParameterValues());
        }

        if (tipData instanceof ContinuousTraitDataModel) {

            setParameter(makeElementaryVector(index), tipData.getParameter());

        } else if (tipData instanceof ElementaryVectorDataModel) {

            ElementaryVectorDataModel elementaryData = (ElementaryVectorDataModel) tipData;

            int tip = index / dimTrait;
            int dim = index % dimTrait;

            elementaryData.setTipTraitDimParameters(tip, 0, dim);

        } else {
            throw new RuntimeException("Not yet implemented");
        }

        double[] column = productProvider.getProduct(parameter);

        if (RESET_DATA) {
            setParameter(savedParameter, parameter);
        }

        return column;
    }

    private WrappedVector makeElementaryVector(int index) {
        double[] vector = new double[numTaxa * dimTrait];
        vector[index] = 1.0;

        return new WrappedVector.Raw(vector);
    }

    private double[] expensiveColumn(int index) {
        return likelihoodDelegate.getTreeTraitPrecision()[index];
    }

    private void debug(double[] result, int index) {
        double[] expensiveResult = expensiveColumn(index);
        System.err.println("via FCD: " + new WrappedVector.Raw(result));
        System.err.println("direct : " + new WrappedVector.Raw(expensiveResult));
        System.err.println();
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        for (int taxon = 0; taxon < numTaxa; ++taxon) {
            for (int trait = 0; trait < dimTrait; ++trait) {
                sb.append(new WrappedVector.Raw(getColumn(taxon * dimTrait + trait))).append("\n");
            }
        }
        return sb.toString();
    }

    private static final boolean DEBUG = false;
    private static final boolean DEBUG_CACHE = false;
    private static final boolean RESET_DATA = false;
}
