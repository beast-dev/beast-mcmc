/*
 * TreeTipGradient.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.ops.CommonOps;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BranchRateGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final int nTraits;
    private final int dim;
    private final Parameter rateParameter;
    private final ArbitraryBranchRates branchRateModel;
    private final MatrixParameterInterface diffusionParameter;

    private final DenseMatrix64F L;
    private final DenseMatrix64F Lt;
    private final DenseMatrix64F matrix0;
    private final DenseMatrix64F matrix1;
    private final DenseMatrix64F vector0;

    public BranchRateGradient(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                              Parameter rateParameter) {

        assert(treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;
        this.diffusionParameter = likelihoodDelegate.getDiffusionModel().getPrecisionParameter();

        BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        this.branchRateModel = (brm instanceof ArbitraryBranchRates) ? (ArbitraryBranchRates) brm : null; 

        // TODO Move into different constructor / parser
        String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        TreeTrait<List<BranchSufficientStatistics>> unchecked = treeDataLikelihood.getTreeTrait(bcdName);
        treeTraitProvider = unchecked;

        assert (treeTraitProvider != null);

        nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }
        dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        assert (dim == diffusionParameter.getColumnDimension());

        if (likelihoodDelegate.getIntegrator() instanceof SafeMultivariateWithDriftIntegrator) {
            throw new RuntimeException("Not yet implemented with drift");
        }

        L = new DenseMatrix64F(dim, dim);
        Lt = new DenseMatrix64F(dim, dim);
        matrix0 = new DenseMatrix64F(dim, dim);
        matrix1 = new DenseMatrix64F(dim, dim);
        vector0 = new DenseMatrix64F(dim, 1);
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        getCholeskyDecomposition(diffusionParameter, L, Lt);

        // TODO Do single call to traitProvider with node == null (get full tree)
//        List<BranchSufficientStatistics> statisticsForTree = (List<BranchSufficientStatistics>)
//                treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {

                List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);

                assert (statisticsForNode.size() == nTraits);

                double gradient = 0.0;

                for (int trait = 0; trait < nTraits; ++trait) {
                    gradient += getGradientForBranch(statisticsForNode.get(trait), L, Lt);
                }

                final int destinationIndex = getParameterIndexFromNode(node);
                assert (destinationIndex != -1);

                result[destinationIndex] = gradient;
            }
        }

        return result;
    }

    private void getCholeskyDecomposition(MatrixParameterInterface input, DenseMatrix64F L, DenseMatrix64F Lt) {

        CholeskyDecomposition<DenseMatrix64F> cholesky = DecompositionFactory.chol(dim, true);
        DenseMatrix64F matrix = DenseMatrix64F.wrap(dim, dim, input.getParameterValues());

        if( !cholesky.decompose(matrix)) {
            throw new RuntimeException("Cholesky decomposition failed!");
        }

//        return cholesky.getT(null);
        cholesky.getT(L);
        CommonOps.transpose(L, Lt);
    }

    private int getParameterIndexFromNode(NodeRef node) {
        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

    private double getGradientForBranch(BranchSufficientStatistics statistics,
                                               DenseMatrix64F L, DenseMatrix64F Lt) {

        final NormalSufficientStatistics child = statistics.getChild();
        final NormalSufficientStatistics branch = statistics.getBranch();
        final NormalSufficientStatistics parent = statistics.getParent();

        if (DEBUG) {
            System.err.println("L = " + NormalSufficientStatistics.toVectorizedString(L));
            System.err.println("B = " + statistics.toVectorizedString());
        }

        // TODO Consider having traitProvider return variance, if rate-limiting
        DenseMatrix64F V = matrix0;
        DenseMatrix64F delta = vector0;

        for (int row = 0; row < dim; ++row) {
            delta.unsafe_set(row, 0,
                    child.getMean(row) - parent.getMean(row) // TODO Correct for drift
            );
        }

        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {

                V.unsafe_set(row, col,
                        child.getVariance(row, col)
                                - 2 * branch.getVariance(row, col)
                                + parent.getVariance(row, col)
                                + delta.unsafe_get(row, 0) * delta.unsafe_get(col, 0)
                );
            }
        }

        if (DEBUG) {
            System.err.println("V = " + NormalSufficientStatistics.toVectorizedString(V));
        }

        DenseMatrix64F tmp = matrix1;
        CommonOps.mult(Lt, V, tmp);
        DenseMatrix64F LtVL = matrix0;
        CommonOps.mult(tmp, L, LtVL);

        double gradient = LtVL.unsafe_get(0, 0);

        if (DEBUG) {
            System.err.println("LtVL = " + NormalSufficientStatistics.toVectorizedString(LtVL));
            System.err.println("gradient = " + gradient);
        }

        return gradient;
    }

    @Override
    public String getReport() {
        return (new dr.math.matrixAlgebra.Vector(getGradientLogDensity())).toString();
    }

    private static final boolean DEBUG = true;
}
