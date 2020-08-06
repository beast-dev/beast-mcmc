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

import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.inference.model.*;
import dr.math.KroneckerOperation;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public class TreeKroneckerPrecisionColumnProvider extends TreePrecisionColumnProvider {

    private final int precisionDim;
    private final MatrixParameterInterface diffusionPrecision;
    private final double[] buffer;
    private final Map<Integer, double[]> kroneckerCache = new HashMap<>();
    private final boolean extendTipBranchTransformed;

    public TreeKroneckerPrecisionColumnProvider(TreePrecisionTraitProductProvider productProvider,
                                                MultivariateDiffusionModel diffusionModel, boolean extendTipBranchTransformed) {

        super(productProvider);

        this.diffusionPrecision = diffusionModel.getPrecisionParameter();
        this.precisionDim = diffusionPrecision.getColumnDimension();
        this.extendTipBranchTransformed = extendTipBranchTransformed;
        addVariable(diffusionPrecision);

        this.buffer = new double[precisionDim];
    }

    @Override
    public double[] getColumn(int index) {

        double[] column = kroneckerCache.get(index);
        if (column == null) {

            if (DEBUG) {
                System.err.println("treeKronecker key " + index + " not found");
            }

            int treeIndex = index / precisionDim;
            int precisionIndex = index % precisionDim;

            double[] treeColumn = super.getColumn(treeIndex);
            double[] precisionColumn = getPrecisionColumn(precisionIndex);

            column = KroneckerOperation.product(
                    treeColumn, 1, treeColumn.length,
                    precisionColumn, 1, precisionColumn.length);

            kroneckerCache.put(index, column);

        } else {
            if (DEBUG) {
                System.err.println("treeKronecker key " + index + " found");
            }
        }
        
        return column;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == diffusionPrecision) {
            kroneckerCache.clear();
        }

        if (extendTipBranchTransformed){
            super.clearTreeCache();
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            kroneckerCache.clear();
        }
        super.handleModelChangedEvent(model, object, index);
    }

    private double[] getPrecisionColumn(int index) {
        double[] column = buffer;

        for (int dim = 0; dim < precisionDim; ++dim) {
            column[dim] = diffusionPrecision.getParameterValue(dim, index);
        }

        return column;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        double[] firstColumn = getColumn(0);
        sb.append(new WrappedVector.Raw(firstColumn)).append("\n");

        for (int index = 1; index < firstColumn.length; ++index) {
            sb.append(new WrappedVector.Raw(getColumn(index))).append("\n");
        }

        return sb.toString();
    }

    private static final boolean DEBUG = false;
}
