/*
 * RepeatedMeasuresTraitDataModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */
public class RepeatedMeasuresTraitDataModel extends ContinuousTraitDataModel implements ContinuousTraitPartialsProvider,
        ModelExtensionProvider.NormalExtensionProvider {

    private final String traitName;
    private final MatrixParameterInterface samplingPrecisionParameter;
    private boolean diagonalOnly = false;
//    private DenseMatrix64F samplingVariance;
    private boolean variableChanged = true;
    private boolean varianceKnown = false;

    private Matrix samplingPrecision;
    private Matrix samplingVariance;
    private Matrix storedSamplingPrecision;
    private Matrix storedSamplingVariance;
    private boolean storedVarianceKnown = false;
    private boolean storedVariableChanged = true;


    public RepeatedMeasuresTraitDataModel(String name,
                                          CompoundParameter parameter,
                                          List<Integer> missingIndices,
//                                          boolean[] missindIndicators,
                                          boolean useMissingIndices,
                                          final int dimTrait,
                                          MatrixParameterInterface samplingPrecision) {
        super(name, parameter, missingIndices, useMissingIndices, dimTrait, PrecisionType.FULL);
        this.traitName = name;
        this.samplingPrecisionParameter = samplingPrecision;
        addVariable(samplingPrecision);

        calculatePrecisionInfo();

//        this.samplingVariance = new Matrix(samplingPrecision.getParameterAsMatrix()).inverse();
        this.samplingVariance = null;


        samplingPrecisionParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                samplingPrecision.getDimension()));

    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        assert (numTraits == 1);
        assert (samplingPrecision.rows() == dimTrait && samplingPrecision.columns() == dimTrait);

        recomputeVariance();

        if (fullyObserved) {
            return new double[dimTrait + 1];
        }

        double[] partial = super.getTipPartial(taxonIndex, fullyObserved);
        DenseMatrix64F V = MissingOps.wrap(partial, dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

        //TODO: remove diagonalOnly part
        if (diagonalOnly) {
            for (int index = 0; index < dimTrait; index++) {
                V.set(index, index, V.get(index, index) + 1 / samplingPrecision.component(index, index));
            }
        } else {
            for (int i = 0; i < dimTrait; i++) {
                for (int j = 0; j < dimTrait; j++) {
                    V.set(i, j, V.get(i, j) + samplingVariance.component(i, j));
                }
            }
        }


        DenseMatrix64F P = new DenseMatrix64F(dimTrait, dimTrait);
        MissingOps.safeInvert2(V, P, false); //TODO this isn't necessary when this is fully observed

        MissingOps.unwrap(P, partial, dimTrait);
        MissingOps.unwrap(V, partial, dimTrait + dimTrait * dimTrait);

        if (DEBUG) {
            System.err.println("taxon " + taxonIndex);
            System.err.println("\tprecision: " + P);
            System.err.println("\tmean: " + new WrappedVector.Raw(partial, 0, dimTrait));
        }

        return partial;
    }


    private void recomputeVariance() {
        checkVariableChanged();
        if (!varianceKnown) {
            samplingVariance = samplingPrecision.inverse();
            varianceKnown = true;
        }
    }

    public Matrix getSamplingVariance() {
        recomputeVariance();
        return samplingVariance;
    }

    public String getTraitName() {
        return traitName;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (variable == samplingPrecisionParameter) {

            variableChanged = true;
            varianceKnown = false;
            fireModelChanged();
        }
    }

    private void calculatePrecisionInfo() {
        samplingPrecision = new Matrix(samplingPrecisionParameter.getParameterAsMatrix());
    }

    private void checkVariableChanged() {
        if (variableChanged) {
            calculatePrecisionInfo();
            variableChanged = false;
            varianceKnown = false;
        }
    }

    @Override
    protected void storeState() {
        storedSamplingPrecision = samplingPrecision.clone();
        storedSamplingVariance = samplingVariance.clone();
        storedVarianceKnown = varianceKnown;
        storedVariableChanged = variableChanged;
    }

    @Override
    protected void restoreState() {
        Matrix tmp = samplingPrecision;
        samplingPrecision = storedSamplingPrecision;
        storedSamplingPrecision = tmp;

        tmp = samplingVariance;
        samplingVariance = storedSamplingVariance;
        storedSamplingVariance = tmp;

        varianceKnown = storedVarianceKnown;
        variableChanged = storedVariableChanged;
    }

    @Override
    public ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate,
                                                            TreeTrait treeTrait, Tree tree) {
        checkVariableChanged();
        return new ContinuousExtensionDelegate.MultivariateNormalExtensionDelegate(delegate, treeTrait,
                this, tree);
    }

    @Override
    public DenseMatrix64F getExtensionVariance() {
        recomputeVariance();
        double[] buffer = samplingVariance.toArrayComponents();
        return DenseMatrix64F.wrap(dimTrait, dimTrait, buffer);
    }

    @Override
    public MatrixParameterInterface getExtensionPrecision() {
        checkVariableChanged();
        return samplingPrecisionParameter;
    }

    @Override
    public double[] transformTreeTraits(double[] treeTraits) {
        return treeTraits;
    }

    @Override
    public int getDataDimension() {
        return dimTrait;
    }

    private static final boolean DEBUG = false;

}
