/*
 * BlombergKStatistic.java
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

package dr.inference.model;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.math.matrixAlgebra.*;

import java.util.Arrays;


public class BlombergKStatistic extends Statistic.Abstract implements ModelListener {

    public static final String BLOMBERGS_K = "blombergsK";

    private final TreeDataLikelihood traitLikelihood;
    private final TreeModel tree;
    private final TreeTrait treeTrait;
    private boolean needToUpdateTree = true;
    private final int traitDim;
    private Matrix Linv;
    private final int treeDim;
    private double expectedRatio;
    private final ContinuousDataLikelihoodDelegate delegate;
    private final double[] k;

    public BlombergKStatistic(TreeDataLikelihood traitLikelihood, String traitName) {
        this.traitLikelihood = traitLikelihood;
        this.tree = (TreeModel) traitLikelihood.getTree();
        tree.addModelListener(this);

        this.treeTrait = traitLikelihood.getTreeTrait(traitName);
        this.traitDim = traitLikelihood.getDataLikelihoodDelegate().getTraitDim();
        this.treeDim = tree.getTaxonCount();
        this.delegate = (ContinuousDataLikelihoodDelegate) traitLikelihood.getDataLikelihoodDelegate();
        this.k = new double[traitDim];
    }


    @Override
    public int getDimension() {
        return traitDim;
    }

    @Override
    public double getStatisticValue(int dim) {
        if (dim == 0) {
            computeStatistics();
        }
        return k[dim];
    }


    public void computeStatistics() {
        if (needToUpdateTree) {
            double[][] treeStructure = MultivariateTraitDebugUtilities.getTreeVariance(tree,
                    traitLikelihood.getBranchRateModel(),
                    1.0, Double.POSITIVE_INFINITY); //TODO: make sure order is right
            //TODO: don't actually need to construct or invert this. can use Ho & Ane 2014 for all calculations


            SymmetricMatrix V = new SymmetricMatrix(treeStructure);
            Matrix L;
            try {
                CholeskyDecomposition chol = new CholeskyDecomposition(V);
                L = new Matrix(chol.getL());
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
                throw new RuntimeException();
            }

            Linv = L.inverse().transpose();
            double[] ones = new double[treeDim];
            Arrays.fill(ones, 1);

            Vector l;

            try {
                l = Linv.product(new Vector(ones));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
                throw new RuntimeException();
            }

            double sumInverse = 0;
            for (int i = 0; i < treeDim; i++) {
                sumInverse += l.component(i) * l.component(i);
            }

            double trace = 0;
            for (int i = 0; i < treeDim; i++) {
                trace += treeStructure[i][i];
            }

            expectedRatio = (trace - treeDim / sumInverse) / (treeDim - 1);
            needToUpdateTree = false;
        }


        double[] treeTraits = (double[]) treeTrait.getTrait(tree, null);


        double[] mean = delegate.getPostOrderRootMean();

        double[] thisTrait = new double[treeDim];

        for (int trait = 0; trait < traitDim; trait++) {
            for (int taxon = 0; taxon < treeDim; taxon++) {
                thisTrait[taxon] = treeTraits[taxon * traitDim + trait] - mean[trait];
            }


            Vector contrasts;

            try {
                contrasts = Linv.product(new Vector(thisTrait));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
                throw new RuntimeException();
            }

            double mse0 = sumSquares(thisTrait);
            double mse = sumSquares(contrasts.toComponents());

            k[trait] = (mse0 / mse) / expectedRatio;
        }
    }

    private double sumSquares(double[] x) {
        double ss = 0;
        for (int i = 0; i < x.length; i++) {
            ss += x[i] * x[i];
        }
        return ss;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        needToUpdateTree = true;
    }

    @Override
    public void modelRestored(Model model) {
        needToUpdateTree = true;
    }
}
