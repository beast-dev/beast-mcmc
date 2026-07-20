package test.dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalAdjointBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalTransitionBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalTransitionAdjointUtils;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Random;

/**
 * Rotated/block-basis OU work, Phase 3 (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md).
 *
 * <p>Compares the analytic gradient of the Phase 2 backend's log-likelihood
 * (computed via {@link BlockBasisCanonicalAdjointBuilder}, working entirely
 * in {@code Ytilde} coordinates plus the explicit Jacobian-term gradient)
 * against a finite-difference gradient of the independent, already-tested
 * natural-basis reference computation (the same one
 * {@code BlockBasisCanonicalTransitionBuilderTest} uses), for the same
 * fixed 3-tip tree {@code ((A:1,B:1):1,C:2)}, fully-observed tips only.
 *
 * <p>Scope: {@code R} (rotation), {@code rho/theta/t/scalar} (the
 * block-diagonal drift's native parameters), and the diffusion matrix
 * {@code Q}. Stationary-mean ({@code beta}) gradients are out of scope —
 * see the plan's Phase 3 scope note.
 */
public class BlockBasisCanonicalAdjointBuilderTest extends TestCase {

    private static final int DIM = 3;
    private static final double FD_DELTA = 1e-6;
    private static final double TOL = 1e-4;
    private static final int TRIALS = 5;
    private static final double KAPPA0 = 10.0;

    public void testAnalyticGradientsMatchFiniteDifference() {
        final Random random = new Random(11);
        for (int trial = 0; trial < TRIALS; ++trial) {
            final MatrixParameter rotationParam = randomWellConditionedRotation(random);
            final Parameter scalarParam = new Parameter.Default(new double[]{0.5 + random.nextDouble()});
            final Parameter rhoParam = new Parameter.Default(new double[]{0.3 + 0.6 * random.nextDouble()});
            final Parameter thetaParam = new Parameter.Default(new double[]{(random.nextDouble() - 0.5) * 1.2});
            final Parameter tParam = new Parameter.Default(new double[]{(random.nextDouble() - 0.5) * 0.4});
            final BlockDiagonalPolarStableMatrixParameter blockParameter =
                    new BlockDiagonalPolarStableMatrixParameter(
                            "block" + trial, rotationParam, scalarParam, rhoParam, thetaParam, tParam);

            final MatrixParameter diffusion = randomSpdMatrix(random);
            final double[] beta = {random.nextGaussian(), random.nextGaussian(), random.nextGaussian()};
            final double[][] tipValues = {
                    {11.0, 12.0, -2.0},
                    {1.0, 2.0, 5.0},
                    {3.0, 3.0, -4.0}
            };

            final AnalyticGradients analytic = computeAnalyticGradients(blockParameter, diffusion, beta, tipValues);

            checkRotationGradient(blockParameter, diffusion, beta, tipValues, rotationParam, analytic.rGrad, trial);
            checkNativeBlockGradient(blockParameter, diffusion, beta, tipValues, scalarParam, analytic.nativeGrad[0], "scalar", trial);
            checkNativeBlockGradient(blockParameter, diffusion, beta, tipValues, rhoParam, analytic.nativeGrad[1], "rho", trial);
            checkNativeBlockGradient(blockParameter, diffusion, beta, tipValues, thetaParam, analytic.nativeGrad[2], "theta", trial);
            checkNativeBlockGradient(blockParameter, diffusion, beta, tipValues, tParam, analytic.nativeGrad[3], "t", trial);
            checkDiffusionGradient(blockParameter, diffusion, beta, tipValues, analytic.qGradFlat, trial);
        }
    }

    private static final class AnalyticGradients {
        final DenseMatrix64F rGrad;
        final double[] nativeGrad;
        final double[] qGradFlat;

        AnalyticGradients(final DenseMatrix64F rGrad, final double[] nativeGrad, final double[] qGradFlat) {
            this.rGrad = rGrad;
            this.nativeGrad = nativeGrad;
            this.qGradFlat = qGradFlat;
        }
    }

    private AnalyticGradients computeAnalyticGradients(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                                        final MatrixParameter diffusion,
                                                        final double[] beta,
                                                        final double[][] tipValues) {
        final BlockBasisGlobalSnapshot snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusion);
        snapshot.ensureCurrent();
        final BlockBasisCanonicalTransitionBuilder builder = new BlockBasisCanonicalTransitionBuilder(
                DIM, snapshot.getBlockStarts(), snapshot.getBlockSizes());
        final BlockBasisCanonicalAdjointBuilder adjointBuilder = new BlockBasisCanonicalAdjointBuilder(
                DIM, snapshot.getBlockStarts(), snapshot.getBlockSizes());

        final CanonicalGaussianTransition transitionLen1 = new CanonicalGaussianTransition(DIM);
        builder.fillCanonicalTransition(snapshot, 1.0, transitionLen1);
        final CanonicalGaussianTransition transitionLen2 = new CanonicalGaussianTransition(DIM);
        builder.fillCanonicalTransition(snapshot, 2.0, transitionLen2);

        final double[][] tipValuesTilde = new double[tipValues.length][DIM];
        for (int t = 0; t < tipValues.length; ++t) {
            transformToTilde(snapshot, tipValues[t], beta, tipValuesTilde[t]);
        }

        final CanonicalGaussianMessageOps.Workspace ws = new CanonicalGaussianMessageOps.Workspace(DIM);

        final CanonicalGaussianState fromA = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValuesTilde[0], fromA);
        final CanonicalGaussianState fromB = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValuesTilde[1], fromB);
        final CanonicalGaussianState n1 = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromA, fromB, n1);

        final CanonicalGaussianState fromN1AtRoot = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.pushBackward(n1, transitionLen1, ws, fromN1AtRoot);
        final CanonicalGaussianState fromC = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen2, tipValuesTilde[2], fromC);
        final CanonicalGaussianState rootMessage = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromN1AtRoot, fromC, rootMessage);

        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();
        final double[] zeroMean = new double[DIM];
        final double[] priorCovFlat = new double[DIM * DIM];
        for (int i = 0; i < DIM * DIM; ++i) {
            priorCovFlat[i] = vTilde.data[i] / KAPPA0;
        }
        final GaussianFormConverter.Workspace converterWs = new GaussianFormConverter.Workspace();
        converterWs.ensureDim(DIM);
        final CanonicalGaussianState priorState = new CanonicalGaussianState(DIM);
        GaussianFormConverter.fillStateFromMoments(zeroMean, priorCovFlat, DIM, converterWs, priorState);

        final CanonicalGaussianState combined = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(rootMessage, priorState, combined);

        // Backward/smoothing pass.
        final CanonicalGaussianState rootExcludingN1 = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromC, priorState, rootExcludingN1);
        final CanonicalGaussianState outsideN1 = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.pushForward(rootExcludingN1, transitionLen1, ws, outsideN1);
        final CanonicalGaussianState n1Posterior = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(n1, outsideN1, n1Posterior);

        final double[] compressedDAcc = new double[blockParameter.getCompressedDDimension()];
        final DenseMatrix64F gSAcc = new DenseMatrix64F(DIM, DIM);
        final DenseMatrix64F rGradAcc = new DenseMatrix64F(DIM, DIM);
        final double[] qGradFlatAcc = new double[DIM * DIM];

        final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWs =
                new CanonicalTransitionAdjointUtils.Workspace(DIM);
        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(DIM);

        // Tips A, B: degenerate-observation Fisher contribution using N1's full marginal posterior.
        accumulateTipPullback(n1Posterior, tipValuesTilde[0], transitionLen1, 1.0,
                transitionAdjointWs, adjoints, adjointBuilder, snapshot, compressedDAcc, gSAcc, rGradAcc);
        accumulateTipPullback(n1Posterior, tipValuesTilde[1], transitionLen1, 1.0,
                transitionAdjointWs, adjoints, adjointBuilder, snapshot, compressedDAcc, gSAcc, rGradAcc);

        // Internal branch N1 -> root: full (non-degenerate) pair posterior.
        final CanonicalGaussianState pairState = new CanonicalGaussianState(2 * DIM);
        CanonicalGaussianMessageOps.buildPairPosterior(rootExcludingN1, transitionLen1, n1, pairState);
        final CanonicalBranchMessageContributionUtils.Workspace pairWs =
                new CanonicalBranchMessageContributionUtils.Workspace(DIM);
        final CanonicalBranchMessageContribution contributionN1Root = new CanonicalBranchMessageContribution(DIM);
        CanonicalBranchMessageContributionUtils.fillFromPairState(pairState, pairWs, contributionN1Root);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transitionLen1, contributionN1Root, transitionAdjointWs, adjoints);
        adjointBuilder.accumulateBranchPullback(snapshot, 1.0, adjoints, compressedDAcc, gSAcc);

        // Tip C: degenerate-observation Fisher contribution using root's full marginal posterior.
        accumulateTipPullback(combined, tipValuesTilde[2], transitionLen2, 2.0,
                transitionAdjointWs, adjoints, adjointBuilder, snapshot, compressedDAcc, gSAcc, rGradAcc);

        // Root prior's direct dependence on Vtilde.
        accumulateRootPriorVtildePullback(combined, priorState, vTilde, adjointBuilder, gSAcc);

        adjointBuilder.finalizeGlobalPullback(snapshot, gSAcc, compressedDAcc, rGradAcc, qGradFlatAcc);
        adjointBuilder.addJacobianTermGradient(snapshot, tipValues.length, rGradAcc);

        final double[] nativeGrad = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedDAcc, nativeGrad);

        return new AnalyticGradients(rGradAcc, nativeGrad, qGradFlatAcc);
    }

    /**
     * Fisher contribution for a tip with degenerate (observed, zero-variance)
     * child value: {@code E[YY^T]=y y^T}, {@code E[XY^T]=E[X] y^T}, etc., with
     * {@code E[X], E[XX^T]} from the parent's full marginal posterior.
     */
    private void accumulateTipPullback(final CanonicalGaussianState parentPosterior,
                                       final double[] yObserved,
                                       final CanonicalGaussianTransition transition,
                                       final double dt,
                                       final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWs,
                                       final CanonicalLocalTransitionAdjoints adjoints,
                                       final BlockBasisCanonicalAdjointBuilder adjointBuilder,
                                       final BlockBasisGlobalSnapshot snapshot,
                                       final double[] compressedDAcc,
                                       final DenseMatrix64F gSAcc,
                                       final DenseMatrix64F rGradAcc) {
        final GaussianFormConverter.Workspace ws = new GaussianFormConverter.Workspace();
        ws.ensureDim(DIM);
        final double[] parentMean = new double[DIM];
        final double[] parentCovFlat = new double[DIM * DIM];
        GaussianFormConverter.fillMomentsFromState(parentPosterior, parentMean, parentCovFlat, DIM, ws);

        final CanonicalBranchMessageContribution contribution = new CanonicalBranchMessageContribution(DIM);
        for (int i = 0; i < DIM; ++i) {
            contribution.dLogL_dInformationX[i] = parentMean[i];
            contribution.dLogL_dInformationY[i] = yObserved[i];
            for (int j = 0; j < DIM; ++j) {
                final double xxMoment = parentCovFlat[i * DIM + j] + parentMean[i] * parentMean[j];
                contribution.dLogL_dPrecisionXX[i * DIM + j] = -0.5 * xxMoment;
                contribution.dLogL_dPrecisionXY[i * DIM + j] = -0.5 * parentMean[i] * yObserved[j];
                contribution.dLogL_dPrecisionYX[i * DIM + j] = -0.5 * yObserved[i] * parentMean[j];
                contribution.dLogL_dPrecisionYY[i * DIM + j] = -0.5 * yObserved[i] * yObserved[j];
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;

        // dLogL/dyObserved, holding the transition fixed: out.precision is y-independent;
        // out.information[i] = informationX[i] - sum_j precisionXY[i,j]*y[j];
        // out.logNormalizer = ... + 0.5*y^T*precisionYY*y - informationY^T*y.
        // dLogL/dy = -precisionXY^T*parentMean - precisionYY*y + informationY (dLogL/d(out.information)=parentMean, dLogL/d(out.logNormalizer)=-1).
        final double[] dLogL_dY = new double[DIM];
        for (int k = 0; k < DIM; ++k) {
            double precisionYYy = 0.0;
            double precisionXYTransposeMean = 0.0;
            for (int i = 0; i < DIM; ++i) {
                precisionYYy += transition.precisionYY[k * DIM + i] * yObserved[i];
                precisionXYTransposeMean += transition.precisionXY[i * DIM + k] * parentMean[i];
            }
            dLogL_dY[k] = -precisionXYTransposeMean - precisionYYy + transition.informationY[k];
        }
        adjointBuilder.addTipTransformRGradient(snapshot, yObserved, dLogL_dY, rGradAcc);

        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition, contribution, transitionAdjointWs, adjoints);
        adjointBuilder.accumulateBranchPullback(snapshot, dt, adjoints, compressedDAcc, gSAcc);
    }

    /** The root prior's own direct dependence on {@code Vtilde} (mean is zero regardless of Vtilde). */
    private void accumulateRootPriorVtildePullback(final CanonicalGaussianState rootPosterior,
                                                    final CanonicalGaussianState priorState,
                                                    final DenseMatrix64F vTilde,
                                                    final BlockBasisCanonicalAdjointBuilder adjointBuilder,
                                                    final DenseMatrix64F gSAcc) {
        final GaussianFormConverter.Workspace ws = new GaussianFormConverter.Workspace();
        ws.ensureDim(DIM);
        final double[] rootMean = new double[DIM];
        final double[] rootCovFlat = new double[DIM * DIM];
        GaussianFormConverter.fillMomentsFromState(rootPosterior, rootMean, rootCovFlat, DIM, ws);

        final double[] dLogL_dP = new double[DIM * DIM];
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                dLogL_dP[i * DIM + j] = -0.5 * (rootCovFlat[i * DIM + j] + rootMean[i] * rootMean[j]);
            }
        }

        final DenseMatrix64F p = new DenseMatrix64F(DIM, DIM, true, priorState.precision);
        final DenseMatrix64F g = new DenseMatrix64F(DIM, DIM, true, dLogL_dP);
        final DenseMatrix64F temp = new DenseMatrix64F(DIM, DIM);
        final DenseMatrix64F pGp = new DenseMatrix64F(DIM, DIM);
        CommonOps.mult(p, g, temp);
        CommonOps.mult(temp, p, pGp);

        final double[] result = new double[DIM * DIM];
        for (int i = 0; i < DIM * DIM; ++i) {
            result[i] = -(1.0 / KAPPA0) * (pGp.data[i] + 0.5 * p.data[i]);
        }
        adjointBuilder.accumulateRootPriorPullback(result, gSAcc);
    }

    private double computeReferenceLogLikelihood(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                                  final MatrixParameter diffusion,
                                                  final double[] beta,
                                                  final double[][] tipValues) {
        final dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization production =
                new dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization(blockParameter);

        final CanonicalGaussianTransition transitionLen1 = new CanonicalGaussianTransition(DIM);
        production.fillCanonicalTransition(diffusion, beta, 1.0, transitionLen1);
        final CanonicalGaussianTransition transitionLen2 = new CanonicalGaussianTransition(DIM);
        production.fillCanonicalTransition(diffusion, beta, 2.0, transitionLen2);

        final CanonicalGaussianMessageOps.Workspace workspace = new CanonicalGaussianMessageOps.Workspace(DIM);

        final CanonicalGaussianState fromA = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValues[0], fromA);
        final CanonicalGaussianState fromB = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValues[1], fromB);

        final CanonicalGaussianState n1 = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromA, fromB, n1);

        final CanonicalGaussianState fromN1AtRoot = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.pushBackward(n1, transitionLen1, workspace, fromN1AtRoot);

        final CanonicalGaussianState fromC = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen2, tipValues[2], fromC);

        final CanonicalGaussianState rootMessage = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromN1AtRoot, fromC, rootMessage);

        final BlockBasisGlobalSnapshot snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusion);
        snapshot.ensureCurrent();
        final DenseMatrix64F vNatural = naturalStationaryCovariance(snapshot);

        final double[] priorCovFlat = new double[DIM * DIM];
        for (int i = 0; i < DIM * DIM; ++i) {
            priorCovFlat[i] = vNatural.data[i] / KAPPA0;
        }
        final GaussianFormConverter.Workspace converterWorkspace = new GaussianFormConverter.Workspace();
        converterWorkspace.ensureDim(DIM);
        final CanonicalGaussianState priorState = new CanonicalGaussianState(DIM);
        GaussianFormConverter.fillStateFromMoments(beta, priorCovFlat, DIM, converterWorkspace, priorState);

        final CanonicalGaussianState combined = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(rootMessage, priorState, combined);
        return CanonicalGaussianMessageOps.normalizationShift(combined, workspace);
    }

    private void checkRotationGradient(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                       final MatrixParameter diffusion,
                                       final double[] beta,
                                       final double[][] tipValues,
                                       final MatrixParameter rotationParam,
                                       final DenseMatrix64F analyticRGrad,
                                       final int trial) {
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                if (i == j) {
                    continue;
                }
                final double original = rotationParam.getParameterValue(i, j);
                rotationParam.setParameterValue(i, j, original + FD_DELTA);
                final double plus = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);
                rotationParam.setParameterValue(i, j, original - FD_DELTA);
                final double minus = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);
                rotationParam.setParameterValue(i, j, original);

                final double numeric = (plus - minus) / (2.0 * FD_DELTA);
                assertEquals("trial " + trial + ": dLogL/dR[" + i + "," + j + "]",
                        numeric, analyticRGrad.get(i, j), TOL);
            }
        }
    }

    private void checkNativeBlockGradient(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                          final MatrixParameter diffusion,
                                          final double[] beta,
                                          final double[][] tipValues,
                                          final Parameter param,
                                          final double analyticGrad,
                                          final String label,
                                          final int trial) {
        final double original = param.getParameterValue(0);
        param.setParameterValue(0, original + FD_DELTA);
        final double plus = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);
        param.setParameterValue(0, original - FD_DELTA);
        final double minus = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);
        param.setParameterValue(0, original);

        final double numeric = (plus - minus) / (2.0 * FD_DELTA);
        assertEquals("trial " + trial + ": dLogL/d" + label, numeric, analyticGrad, TOL);
    }

    private void checkDiffusionGradient(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                        final MatrixParameter diffusion,
                                        final double[] beta,
                                        final double[][] tipValues,
                                        final double[] analyticQGradFlat,
                                        final int trial) {
        for (int i = 0; i < DIM; ++i) {
            for (int j = i; j < DIM; ++j) {
                final double originalIJ = diffusion.getParameterValue(i, j);
                final double originalJI = diffusion.getParameterValue(j, i);

                diffusion.setParameterValue(i, j, originalIJ + FD_DELTA);
                diffusion.setParameterValue(j, i, originalJI + FD_DELTA);
                final double plus = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);

                diffusion.setParameterValue(i, j, originalIJ - FD_DELTA);
                diffusion.setParameterValue(j, i, originalJI - FD_DELTA);
                final double minus = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);

                diffusion.setParameterValue(i, j, originalIJ);
                diffusion.setParameterValue(j, i, originalJI);

                final double numeric = (plus - minus) / (2.0 * FD_DELTA);
                final double analytic = (i == j)
                        ? analyticQGradFlat[i * DIM + j]
                        : analyticQGradFlat[i * DIM + j] + analyticQGradFlat[j * DIM + i];
                assertEquals("trial " + trial + ": dLogL/dQ[" + i + "," + j + "]", numeric, analytic, TOL);
            }
        }
    }

    private static void transformToTilde(final BlockBasisGlobalSnapshot snapshot,
                                         final double[] y,
                                         final double[] beta,
                                         final double[] out) {
        final DenseMatrix64F t = snapshot.getRotationInverse();
        for (int i = 0; i < DIM; ++i) {
            double sum = 0.0;
            for (int j = 0; j < DIM; ++j) {
                sum += t.get(i, j) * (y[j] - beta[j]);
            }
            out[i] = sum;
        }
    }

    private static DenseMatrix64F naturalStationaryCovariance(final BlockBasisGlobalSnapshot snapshot) {
        final DenseMatrix64F r = snapshot.getRotation();
        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();
        final DenseMatrix64F temp = new DenseMatrix64F(DIM, DIM);
        final DenseMatrix64F out = new DenseMatrix64F(DIM, DIM);
        CommonOps.mult(r, vTilde, temp);
        CommonOps.multTransB(temp, r, out);
        return out;
    }

    private static MatrixParameter randomWellConditionedRotation(final Random random) {
        final MatrixParameter rotation = new MatrixParameter("R", DIM, DIM);
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                final double offDiagonal = (i == j) ? 0.0 : (random.nextDouble() - 0.5) * 0.4;
                rotation.setParameterValue(i, j, (i == j ? 1.0 : 0.0) + offDiagonal);
            }
        }
        return rotation;
    }

    private static MatrixParameter randomSpdMatrix(final Random random) {
        final double[][] a = new double[DIM][DIM];
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                a[i][j] = random.nextDouble() - 0.5;
            }
        }
        final MatrixParameter q = new MatrixParameter("Q", DIM, DIM);
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                double sum = 0.0;
                for (int k = 0; k < DIM; ++k) {
                    sum += a[i][k] * a[j][k];
                }
                if (i == j) {
                    sum += 0.5;
                }
                q.setParameterValue(i, j, sum);
            }
        }
        return q;
    }
}
