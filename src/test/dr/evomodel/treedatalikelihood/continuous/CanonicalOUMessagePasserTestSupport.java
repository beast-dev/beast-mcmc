/*
 * CanonicalOUMessagePasserTestSupport.java
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

package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalConjugateRootPriorAdapter;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.HomogeneousCanonicalOUBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.MatrixUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.SequentialCanonicalOUMessagePasser;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

abstract class CanonicalOUMessagePasserTestSupport extends ContinuousTraitTest {

    protected static final double FD_DELTA = 1e-5;
    protected static final double TOL = 1e-4;
    protected static final double BRANCH_TOL = 1e-2;

    protected static final double[][] A_VALUES = {
            {2.0, 0.5, 0.1},
            {0.5, 3.0, 0.3},
            {0.1, 0.3, 1.5}
    };

    protected static final double[] MU_VALUES = {1.0, -0.5, 2.0};
    protected static final double ORTHOGONAL_SCALAR = 1.3;
    protected static final double ORTHOGONAL_RHO = 0.85;
    protected static final double ORTHOGONAL_THETA = 0.25;
    protected static final double ORTHOGONAL_T = -0.08;
    protected static final double[] ORTHOGONAL_ANGLES = {0.2, -0.15, 0.1};

    protected CanonicalOUMessagePasserTestSupport(final String name) {
        super(name);
    }

    @Override
    public void setUp() throws Exception {
        super.setUp();
    }

    protected void checkGradientQ(final OUSetup setup, final String label) {
        final int d2 = dimTrait * dimTrait;

        refreshMessages(setup);
        final double[] gradAnalytic = new double[d2];
        computeJointGradientQ(setup, gradAnalytic, d2);

        final double[] gradNumeric = new double[d2];
        final double[] gradAnalyticSymmetric = new double[d2];
        for (int row = 0; row < dimTrait; row++) {
            for (int col = row; col < dimTrait; col++) {
                final int ij = row * dimTrait + col;
                final int ji = col * dimTrait + row;
                final double numeric;

                if (row == col) {
                    final double original = setup.qMatrix.getParameterValue(row, col);

                    setup.qMatrix.setParameterValue(row, col, original + FD_DELTA);
                    syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
                    final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

                    setup.qMatrix.setParameterValue(row, col, original - FD_DELTA);
                    syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
                    final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

                    setup.qMatrix.setParameterValue(row, col, original);
                    syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
                    numeric = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
                    gradAnalyticSymmetric[ij] = gradAnalytic[ij];
                } else {
                    final double originalIJ = setup.qMatrix.getParameterValue(row, col);
                    final double originalJI = setup.qMatrix.getParameterValue(col, row);

                    setup.qMatrix.setParameterValue(row, col, originalIJ + FD_DELTA);
                    setup.qMatrix.setParameterValue(col, row, originalJI + FD_DELTA);
                    syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
                    final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

                    setup.qMatrix.setParameterValue(row, col, originalIJ - FD_DELTA);
                    setup.qMatrix.setParameterValue(col, row, originalJI - FD_DELTA);
                    syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
                    final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

                    setup.qMatrix.setParameterValue(row, col, originalIJ);
                    setup.qMatrix.setParameterValue(col, row, originalJI);
                    syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
                    numeric = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
                    gradAnalyticSymmetric[ij] = gradAnalytic[ij] + gradAnalytic[ji];
                    gradAnalyticSymmetric[ji] = gradAnalyticSymmetric[ij];
                }

                gradNumeric[ij] = numeric;
                gradNumeric[ji] = numeric;
            }
        }

        System.out.printf("  Analytic sym ∂logL/∂Q: %s%n", Arrays.toString(gradAnalyticSymmetric));
        System.out.printf("  Numeric  sym ∂logL/∂Q: %s%n", Arrays.toString(gradNumeric));
        for (int row = 0; row < dimTrait; row++) {
            for (int col = row; col < dimTrait; col++) {
                final int ij = row * dimTrait + col;
                assertEquals("∂logL/∂Qsym[" + row + "," + col + "] (" + label + ")",
                        gradNumeric[ij], gradAnalyticSymmetric[ij], TOL);
            }
        }
    }

    protected void checkGradientA(final OUSetup setup, final String label) {
        final int d2 = dimTrait * dimTrait;

        refreshMessages(setup);
        final double[] gradAnalytic = new double[d2];
        computeJointGradientA(setup, gradAnalytic);

        final double[] gradNumeric = new double[d2];
        for (int ij = 0; ij < d2; ij++) {
            final int row = ij / dimTrait;
            final int col = ij % dimTrait;
            final double original = setup.aMatrix.getParameterValue(row, col);

            setup.aMatrix.setParameterValue(row, col, original + FD_DELTA);
            final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            setup.aMatrix.setParameterValue(row, col, original - FD_DELTA);
            final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            setup.aMatrix.setParameterValue(row, col, original);
            gradNumeric[ij] = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        }

        System.out.printf("  Analytic ∂logL/∂A: %s%n", Arrays.toString(gradAnalytic));
        System.out.printf("  Numeric  ∂logL/∂A: %s%n", Arrays.toString(gradNumeric));
        for (int ij = 0; ij < d2; ij++) {
            assertEquals("∂logL/∂A[" + (ij / dimTrait) + "," + (ij % dimTrait) + "] (" + label + ")",
                    gradNumeric[ij], gradAnalytic[ij], TOL);
        }
    }

    protected void checkGradientMu(final OUSetup setup, final String label) {
        refreshMessages(setup);
        final double[] gradAnalytic = new double[dimTrait];
        computeJointGradientMu(setup, gradAnalytic, dimTrait * dimTrait);

        final double[] gradNumeric = new double[dimTrait];
        for (int i = 0; i < dimTrait; i++) {
            final double original = setup.muParam.getParameterValue(i);

            setup.muParam.setParameterValue(i, original + FD_DELTA);
            final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            setup.muParam.setParameterValue(i, original - FD_DELTA);
            final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            setup.muParam.setParameterValue(i, original);
            gradNumeric[i] = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        }

        System.out.printf("  Analytic ∂logL/∂mu: %s%n", Arrays.toString(gradAnalytic));
        System.out.printf("  Numeric  ∂logL/∂mu: %s%n", Arrays.toString(gradNumeric));
        for (int i = 0; i < dimTrait; i++) {
            assertEquals("∂logL/∂mu[" + i + "] (" + label + ")",
                    gradNumeric[i], gradAnalytic[i], TOL);
        }
    }

    protected void checkGradientBranchLengths(final OUSetup setup, final String label) {
        refreshMessages(setup);
        final int nodeCount = treeModel.getNodeCount();
        final int rootIndex = treeModel.getRoot().getNumber();

        final double[] gradAnalytic = new double[nodeCount];
        setup.passer.computeGradientBranchLengths(setup.provider, gradAnalytic);

        final double[] gradNumeric = new double[nodeCount];
        for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            if (nodeIndex == rootIndex) {
                continue;
            }

            final NodeRef node = treeModel.getNode(nodeIndex);
            final double baseBranchLength = treeModel.getBranchLength(node);

            setup.branchScale[nodeIndex] += FD_DELTA;
            final double logLPlus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            setup.branchScale[nodeIndex] -= 2.0 * FD_DELTA;
            final double logLMinus = setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);

            setup.branchScale[nodeIndex] += FD_DELTA;
            gradNumeric[nodeIndex] = (logLPlus - logLMinus) / (2.0 * FD_DELTA * baseBranchLength);
        }

        System.out.printf("  Analytic ∂logL/∂t: %s%n", Arrays.toString(gradAnalytic));
        System.out.printf("  Numeric  ∂logL/∂t: %s%n", Arrays.toString(gradNumeric));
        for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            if (nodeIndex == rootIndex) {
                continue;
            }
            assertEquals("∂logL/∂t[" + nodeIndex + "] (" + label + ")",
                    gradNumeric[nodeIndex], gradAnalytic[nodeIndex], BRANCH_TOL);
        }
    }

    protected void refreshMessages(final OUSetup setup) {
        syncPrecisionFromVariance(setup.qMatrix, setup.precisionMatrix);
        setup.passer.computePostOrderLogLikelihood(setup.provider, setup.rootPrior);
        setup.passer.computePreOrder(setup.provider, setup.rootPrior);
    }

    protected void computeJointGradientA(final OUSetup setup, final double[] gradA) {
        setup.passer.computeJointGradients(
                setup.provider,
                gradA,
                new double[dimTrait * dimTrait],
                new double[dimTrait]);
    }

    protected void computeJointGradientQ(final OUSetup setup,
                                         final double[] gradQ,
                                         final int selectionGradientDimension) {
        setup.passer.computeJointGradients(
                setup.provider,
                new double[selectionGradientDimension],
                gradQ,
                new double[dimTrait]);
    }

    protected void computeJointGradientMu(final OUSetup setup,
                                          final double[] gradMu,
                                          final int selectionGradientDimension) {
        setup.passer.computeJointGradients(
                setup.provider,
                new double[selectionGradientDimension],
                new double[dimTrait * dimTrait],
                gradMu);
    }

    protected void checkOrthogonalBlockLikelihoodMatchesDense(final OrthogonalBlockOUSetup setup,
                                                              final String label) {
        final double orthogonalLogLikelihood =
                setup.orthogonal.passer.computePostOrderLogLikelihood(setup.orthogonal.provider, setup.orthogonal.rootPrior);
        final double denseLogLikelihood =
                setup.dense.passer.computePostOrderLogLikelihood(setup.dense.provider, setup.dense.rootPrior);
        System.out.printf("  Orthogonal logL: %.12f%n", orthogonalLogLikelihood);
        System.out.printf("  Dense      logL: %.12f%n", denseLogLikelihood);
        assertEquals("orthogonal-block vs dense log-likelihood (" + label + ")",
                denseLogLikelihood, orthogonalLogLikelihood, 1e-8);
    }

    protected void checkOrthogonalBlockNativeSelectionGradient(final OrthogonalBlockOUSetup setup,
                                                               final String label) {
        refreshMessages(setup.orthogonal);
        final int dimension = setup.nativeParameter.getDimension();
        final double[] gradAnalytic = new double[dimension];
        computeJointGradientA(setup.orthogonal, gradAnalytic);

        final double[] gradNumeric = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            final double original = setup.nativeParameter.getParameterValue(i);

            setup.nativeParameter.setParameterValue(i, original + FD_DELTA);
            final double logLPlus =
                    setup.orthogonal.passer.computePostOrderLogLikelihood(setup.orthogonal.provider, setup.orthogonal.rootPrior);

            setup.nativeParameter.setParameterValue(i, original - FD_DELTA);
            final double logLMinus =
                    setup.orthogonal.passer.computePostOrderLogLikelihood(setup.orthogonal.provider, setup.orthogonal.rootPrior);

            setup.nativeParameter.setParameterValue(i, original);
            gradNumeric[i] = (logLPlus - logLMinus) / (2.0 * FD_DELTA);
        }

        System.out.printf("  Analytic native ∂logL/∂A: %s%n", Arrays.toString(gradAnalytic));
        System.out.printf("  Numeric  native ∂logL/∂A: %s%n", Arrays.toString(gradNumeric));
        for (int i = 0; i < dimension; i++) {
            assertEquals("native orthogonal-block ∂logL/∂A[" + i + "] (" + label + ")",
                    gradNumeric[i], gradAnalytic[i], TOL);
        }
    }

    protected void checkOrthogonalBlockGradientQMatchesDense(final OrthogonalBlockOUSetup setup,
                                                             final String label) {
        refreshMessages(setup.orthogonal);
        refreshMessages(setup.dense);

        final int d2 = dimTrait * dimTrait;
        final double[] orthogonalGradient = new double[d2];
        final double[] denseGradient = new double[d2];
        computeJointGradientQ(setup.orthogonal, orthogonalGradient, setup.nativeParameter.getDimension());
        computeJointGradientQ(setup.dense, denseGradient, d2);

        System.out.printf("  Orthogonal ∂logL/∂Q: %s%n", Arrays.toString(orthogonalGradient));
        System.out.printf("  Dense      ∂logL/∂Q: %s%n", Arrays.toString(denseGradient));
        for (int i = 0; i < d2; i++) {
            assertEquals("orthogonal-block vs dense ∂logL/∂Q[" + i + "] (" + label + ")",
                    denseGradient[i], orthogonalGradient[i], TOL);
        }
    }

    protected OUSetup withPasserParallelism(final OUSetup baseSetup,
                                            final CanonicalTipObservation[] tips,
                                            final int branchGradientParallelism) {
        final SequentialCanonicalOUMessagePasser passer =
                new SequentialCanonicalOUMessagePasser(treeModel, dimTrait, branchGradientParallelism);
        for (int tipIndex = 0; tipIndex < tips.length; tipIndex++) {
            passer.setTipObservation(tipIndex, tips[tipIndex]);
        }
        return new OUSetup(
                baseSetup.aMatrix,
                baseSetup.qMatrix,
                baseSetup.precisionMatrix,
                baseSetup.muParam,
                baseSetup.branchScale,
                baseSetup.provider,
                passer,
                baseSetup.rootPrior);
    }

    protected OUSetup buildOUSetup(final String tag, final CanonicalTipObservation[] tips) {
        final int d = dimTrait;

        final MatrixParameter aMatrix = new MatrixParameter("A." + tag, d, d);
        fillMatrixParameter(aMatrix, A_VALUES);

        return buildOUSetup(tag, tips, aMatrix, MU_VALUES.clone(), this.rootPrior);
    }

    protected OUSetup buildOUSetup(final String tag,
                                   final CanonicalTipObservation[] tips,
                                   final dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior rootPrior) {
        final int d = dimTrait;
        final MatrixParameter aMatrix = new MatrixParameter("A." + tag, d, d);
        fillMatrixParameter(aMatrix, A_VALUES);

        return buildOUSetup(tag, tips, aMatrix, MU_VALUES.clone(), rootPrior);
    }

    protected OUSetup buildOUSetup(final String tag,
                                   final CanonicalTipObservation[] tips,
                                   final MatrixParameter aMatrix,
                                   final double[] muValues) {
        return buildOUSetup(tag, tips, aMatrix, muValues, this.rootPrior);
    }

    protected OUSetup buildOUSetup(final String tag,
                                   final CanonicalTipObservation[] tips,
                                   final MatrixParameter aMatrix,
                                   final double[] muValues,
                                   final dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior rootPrior) {
        final MultivariateElasticModel elasticModel = new MultivariateElasticModel(aMatrix);
        return buildOUSetupFromElasticModel(tag, tips, elasticModel, aMatrix, muValues, rootPrior);
    }

    protected OUSetup buildOUSetupFromElasticModel(final String tag,
                                                   final CanonicalTipObservation[] tips,
                                                   final MultivariateElasticModel elasticModel,
                                                   final MatrixParameter aMatrix,
                                                   final double[] muValues) {
        return buildOUSetupFromElasticModel(tag, tips, elasticModel, aMatrix, muValues, this.rootPrior);
    }

    protected OUSetup buildOUSetupFromElasticModel(final String tag,
                                                   final CanonicalTipObservation[] tips,
                                                   final MultivariateElasticModel elasticModel,
                                                   final MatrixParameter aMatrix,
                                                   final double[] muValues,
                                                   final dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior rootPrior) {
        final int d = dimTrait;
        final int d2 = d * d;

        final double[][] basePrecision = diffusionModel.getPrecisionMatrix();
        final double[] precisionFlat = new double[d2];
        for (int i = 0; i < d; i++) {
            System.arraycopy(basePrecision[i], 0, precisionFlat, i * d, d);
        }
        final double[] qFlat = new double[d2];
        MatrixUtils.invertSymmetric(precisionFlat, qFlat, d);

        final MatrixParameter qMatrix = new MatrixParameter("Q." + tag, d, d);
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                qMatrix.setParameterValueQuietly(i, j, qFlat[i * d + j]);
            }
        }
        qMatrix.fireParameterChangedEvent();

        final MatrixParameter precisionMatrix = new MatrixParameter("precision." + tag, d, d);
        syncPrecisionFromVariance(qMatrix, precisionMatrix);

        final Parameter muParam = new Parameter.Default("mu." + tag, muValues);
        final MultivariateDiffusionModel diffusionModel = new MultivariateDiffusionModel(precisionMatrix);

        final double[] branchScale = new double[treeModel.getNodeCount()];
        Arrays.fill(branchScale, 1.0);
        final ScaledBranchRateModel scaledRateModel = new ScaledBranchRateModel(branchScale);

        final HomogeneousCanonicalOUBranchTransitionProvider provider =
                new HomogeneousCanonicalOUBranchTransitionProvider(
                        treeModel, elasticModel, diffusionModel, muParam, scaledRateModel);
        final SequentialCanonicalOUMessagePasser passer =
                new SequentialCanonicalOUMessagePasser(treeModel, dimTrait);

        for (int tipIndex = 0; tipIndex < tips.length; tipIndex++) {
            passer.setTipObservation(tipIndex, tips[tipIndex]);
        }

        final CanonicalRootPrior canonicalRootPrior =
                new CanonicalConjugateRootPriorAdapter(rootPrior, dimTrait);
        return new OUSetup(aMatrix, qMatrix, precisionMatrix, muParam, branchScale, provider, passer, canonicalRootPrior);
    }

    protected OrthogonalBlockOUSetup buildOrthogonalBlockOUSetup(final String tag,
                                                                 final CanonicalTipObservation[] tips) {
        final Parameter angles = new Parameter.Default("angles." + tag, ORTHOGONAL_ANGLES.clone());
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("R." + tag, dimTrait, angles);
        final Parameter scalar = new Parameter.Default("scalar." + tag, new double[]{ORTHOGONAL_SCALAR});
        final Parameter rho = new Parameter.Default("rho." + tag, new double[]{ORTHOGONAL_RHO});
        final Parameter theta = new Parameter.Default("theta." + tag, new double[]{ORTHOGONAL_THETA});
        final Parameter t = new Parameter.Default("t." + tag, new double[]{ORTHOGONAL_T});

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "A.block." + tag, rotation, scalar, rho, theta, t);

        final OUSetup orthogonalSetup = buildOUSetupFromElasticModel(
                "orth." + tag,
                tips,
                new MultivariateElasticModel(blockSelection),
                null,
                MU_VALUES.clone());

        final MatrixParameter denseSelection = new MatrixParameter("A.dense." + tag, dimTrait, dimTrait);
        fillMatrixParameter(denseSelection, blockSelection.getParameterAsMatrix());
        final OUSetup denseSetup = buildOUSetup(
                "dense." + tag,
                tips,
                denseSelection,
                MU_VALUES.clone());

        return new OrthogonalBlockOUSetup(blockSelection, blockSelection.getParameter(), orthogonalSetup, denseSetup);
    }

    protected CanonicalTipObservation[] buildFullyObservedTips() {
        final Map<String, double[]> traitValues = buildStandardTraitValues();
        final CanonicalTipObservation[] tips = new CanonicalTipObservation[treeModel.getExternalNodeCount()];

        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            final NodeRef node = treeModel.getExternalNode(i);
            final int nodeIndex = node.getNumber();
            if (nodeIndex >= tips.length) {
                throw new IllegalStateException("External node indices must occupy [0, tipCount). Found " + nodeIndex);
            }

            final String taxonId = treeModel.getNodeTaxon(node).getId();
            final double[] values = traitValues.get(taxonId);
            if (values == null) {
                throw new IllegalArgumentException("Missing canonical tip data for taxon " + taxonId);
            }

            final CanonicalTipObservation observation = new CanonicalTipObservation(dimTrait);
            observation.setObserved(values);
            tips[nodeIndex] = observation;
        }

        return tips;
    }

    protected CanonicalTipObservation[] buildPartiallyObservedTips() {
        final CanonicalTipObservation[] tips = buildFullyObservedTips();

        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            final NodeRef node = treeModel.getExternalNode(i);
            final int nodeIndex = node.getNumber();
            final String taxonId = treeModel.getNodeTaxon(node).getId();

            if ("chimp".equals(taxonId)) {
                tips[nodeIndex].setMissing();
                continue;
            }

            if ("bonobo".equals(taxonId)) {
                final boolean[] observedMask = new boolean[]{true, false, true};
                tips[nodeIndex].setPartiallyObserved(tips[nodeIndex].values, observedMask);
            }
        }

        return tips;
    }

    protected Map<String, double[]> buildStandardTraitValues() {
        final Map<String, double[]> out = new HashMap<String, double[]>();
        out.put("human", new double[]{-1.0, 2.0, 3.0});
        out.put("chimp", new double[]{10.0, 12.0, 14.0});
        out.put("bonobo", new double[]{0.5, -2.0, 5.5});
        out.put("gorilla", new double[]{2.0, 5.0, -8.0});
        out.put("orangutan", new double[]{11.0, 1.0, -1.5});
        out.put("siamang", new double[]{1.0, 2.5, 4.0});
        return out;
    }

    protected static void fillMatrixParameter(final MatrixParameter matrix, final double[][] values) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                matrix.setParameterValueQuietly(i, j, values[i][j]);
            }
        }
        matrix.fireParameterChangedEvent();
    }

    protected static void syncPrecisionFromVariance(final MatrixParameter qMatrix,
                                                    final MatrixParameter precisionMatrix) {
        final int dim = qMatrix.getRowDimension();
        final DenseMatrix64F variance = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F precision = new DenseMatrix64F(dim, dim);

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                variance.set(i, j, qMatrix.getParameterValue(i, j));
            }
        }

        if (!CommonOps.invert(variance, precision)) {
            throw new IllegalArgumentException("Variance matrix is singular during numerical gradient check.");
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                precisionMatrix.setParameterValueQuietly(i, j, precision.get(i, j));
            }
        }
        precisionMatrix.fireParameterChangedEvent();
    }

    protected static final class OUSetup {
        final MatrixParameter aMatrix;
        final MatrixParameter qMatrix;
        final MatrixParameter precisionMatrix;
        final Parameter muParam;
        final double[] branchScale;
        final HomogeneousCanonicalOUBranchTransitionProvider provider;
        final SequentialCanonicalOUMessagePasser passer;
        final CanonicalRootPrior rootPrior;

        OUSetup(final MatrixParameter aMatrix,
                final MatrixParameter qMatrix,
                final MatrixParameter precisionMatrix,
                final Parameter muParam,
                final double[] branchScale,
                final HomogeneousCanonicalOUBranchTransitionProvider provider,
                final SequentialCanonicalOUMessagePasser passer,
                final CanonicalRootPrior rootPrior) {
            this.aMatrix = aMatrix;
            this.qMatrix = qMatrix;
            this.precisionMatrix = precisionMatrix;
            this.muParam = muParam;
            this.branchScale = branchScale;
            this.provider = provider;
            this.passer = passer;
            this.rootPrior = rootPrior;
        }
    }

    protected static final class OrthogonalBlockOUSetup {
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection;
        final Parameter nativeParameter;
        final OUSetup orthogonal;
        final OUSetup dense;

        private OrthogonalBlockOUSetup(final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection,
                                       final Parameter nativeParameter,
                                       final OUSetup orthogonal,
                                       final OUSetup dense) {
            this.blockSelection = blockSelection;
            this.nativeParameter = nativeParameter;
            this.orthogonal = orthogonal;
            this.dense = dense;
        }
    }

    protected static final class ScaledBranchRateModel extends AbstractBranchRateModel {

        private final double[] branchScale;

        private ScaledBranchRateModel(final double[] branchScale) {
            super("scaledCanonicalOuTestRateModel");
            this.branchScale = branchScale;
        }

        @Override
        public double getBranchRate(final Tree tree, final NodeRef node) {
            return branchScale[node.getNumber()];
        }

        @Override
        protected void handleModelChangedEvent(final Model model, final Object object, final int index) { }

        @Override
        protected void handleVariableChangedEvent(final Variable variable,
                                                  final int index,
                                                  final Parameter.ChangeType type) { }

        @Override
        protected void storeState() { }

        @Override
        protected void restoreState() { }

        @Override
        protected void acceptState() { }
    }
}
