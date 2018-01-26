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
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.NewTipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import java.util.List;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
public class TreePrecisionTraitProductProvider implements PrecisionMatrixVectorProductProvider, Reportable {

    private final TreeTrait<List<NormalSufficientStatistics>> fullConditionalDensity;
    private final Tree tree;
    private final Parameter dataParameter;

    private final static boolean DEBUG = true;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    public TreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                             ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                             String traitName) {

        String fcdName = NewTipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addNewFullConditionalDensityTrait(traitName);
        }

        this.fullConditionalDensity = castTreeTrait(treeDataLikelihood.getTreeTrait(fcdName));
        this.tree = treeDataLikelihood.getTree();
        this.dataParameter = likelihoodDelegate.getDataModel().getParameter();
        this.likelihoodDelegate = likelihoodDelegate;
    }

    @Override
    public double[] getProduct(Parameter vector) {

        if (vector != dataParameter) {
            throw new IllegalArgumentException("May only compute for trait data vector");
        }

        double[] result = new double[vector.getDimension()];

        List<NormalSufficientStatistics> statistics = fullConditionalDensity.getTrait(tree, null);
        statistics.size();

        // TODO Fill in values in result

        if (DEBUG) {

            double[][] treeTraitPrecision = likelihoodDelegate.getTreeTraitPrecision();
            double[] result2 = new double[result.length];
            for (int row = 0; row < result.length; ++row) {
                double sum = 0.0;
                for (int col = 0; col < treeTraitPrecision[row].length; ++col) {
                    sum += treeTraitPrecision[row][col] * vector.getParameterValue(col);
                }
                result2[row] = sum;
            }

            System.err.println("via FCD: " + new WrappedVector.Raw(result));
            System.err.println("direct : " + new WrappedVector.Raw(result2));
            System.err.println();

            result = result2;
        }

        return result;
    }

    @Override
    public String getReport() {
        double[] result = getProduct(dataParameter);
        return (new WrappedVector.Raw(result, 0, result.length)).toString();
    }

    @SuppressWarnings("unchecked")
    private TreeTrait<List<NormalSufficientStatistics>> castTreeTrait(TreeTrait trait) {
        return trait;
    }
}
