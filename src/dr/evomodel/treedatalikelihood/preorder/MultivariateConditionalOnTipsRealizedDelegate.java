package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 */
public class MultivariateConditionalOnTipsRealizedDelegate extends ConditionalOnTipsRealizedDelegate {

    private static final boolean DEBUG = false;

    final private PartiallyMissingInformation missingInformation;

    public MultivariateConditionalOnTipsRealizedDelegate(String name, Tree tree,
                                                         MultivariateDiffusionModel diffusionModel,
                                                         ContinuousTraitPartialsProvider dataModel,
                                                         ConjugateRootTraitPrior rootPrior,
                                                         ContinuousRateTransformation rateTransformation,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        missingInformation = new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate);
    }

    @Override
    protected void simulateTraitForRoot(final int offsetSample, final int offsetPartial) {

        // TODO Bad programming -- should not need to know about internal layout
        // Layout, offset, dim
        // trait, 0, dT
        // precision, dT, dT * dT
        // variance, dT + dT * dT, dT * dT
        // scalar, dT + 2 * dT * dT, 1

        // Integrate out against prior
        final DenseMatrix64F rootPrec = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);
        final DenseMatrix64F priorPrec = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.mult(Pd, wrap(partialPriorBuffer, offsetPartial + dimTrait, dimTrait, dimTrait), priorPrec);

        final DenseMatrix64F totalPrec = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.add(rootPrec, priorPrec, totalPrec);

        final DenseMatrix64F totalVar = new DenseMatrix64F(dimTrait, dimTrait);
        safeInvert(totalPrec, totalVar, false);

        final double[] tmp = new double[dimTrait];
        final double[] mean = new double[dimTrait];

        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += rootPrec.unsafe_get(g, h) * partialNodeBuffer[offsetPartial + h];
                sum += priorPrec.unsafe_get(g, h) * partialPriorBuffer[offsetPartial + h];
            }
            tmp[g] = sum;
        }
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += totalVar.unsafe_get(g, h) * tmp[h];
            }
            mean[g] = sum;
        }

        final double[][] cholesky = getCholeskyOfVariance(totalVar.getData(), dimTrait);

        MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                mean, 0, // input mean
                cholesky, 1.0, // input variance
                sample, offsetSample, // output sample
                tmpEpsilon);

        if (DEBUG) {
            System.err.println("Attempt to simulate root");
            System.err.println("Root  mean: " + new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait));
            System.err.println("Root  prec: " + rootPrec);
            System.err.println("Prior mean: " + new WrappedVector.Raw(partialPriorBuffer, offsetPartial, dimTrait));
            System.err.println("Prior prec: " + priorPrec);
            System.err.println("Total prec: " + totalPrec);
            System.err.println("Total  var: " + totalVar);


            System.err.println("draw mean: " + new WrappedVector.Raw(mean, 0, dimTrait));
            System.err.println("sample: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));
        }
    }

//    boolean extremeValue(final DenseMatrix64F x) {
//        return extremeValue(new WrappedVector.Raw(x.getData(), 0, x.getNumElements()));
//    }
//
//    boolean extremeValue(final WrappedVector x) {
//        boolean valid = true;
//        for (int i = 0; i < x.getDim() && (valid); ++i) {
//            if (Double.isNaN(x.get(i)) || Math.abs(x.get(i)) > 1E2) {
//                valid = false;
//            }
//        }
//        return !valid;
//    }

    @Override
    protected void simulateTraitForNode(final int nodeIndex,
                                        final int traitIndex,
                                        final int offsetSample,
                                        final int offsetParent,
                                        final int offsetPartial,
                                        final int isExternal,
                                        final double branchPrecision) {

        if (isExternal == 1) { // Is external
            simulateTraitForExternalNode(nodeIndex, traitIndex, offsetSample, offsetParent, offsetPartial, branchPrecision);
        } else {
            simulateTraitForInternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);
        }
    }

    private void  simulateTraitForExternalNode(final int nodeIndex,
                                              final int traitIndex,
                                              final int offsetSample,
                                              final int offsetParent,
                                              final int offsetPartial,
                                              final double branchPrecision) {

        final DenseMatrix64F P0 = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);
        final int missingCount = countFiniteDiagonals(P0);

        if (missingCount == 0) { // Completely observed

            System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);

        } else {

            final int zeroCount = countZeroDiagonals(P0);
            if (zeroCount == dimTrait) { //  All missing completely at random

                final double sqrtScale = Math.sqrt(1.0 / branchPrecision);

                // TODO Drift?
                assert (!likelihoodDelegate.getDiffusionProcessDelegate().hasDrift());

                MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                        sample, offsetParent, // input mean
                        cholesky, sqrtScale, // input variance
                        sample, offsetSample, // output sample
                        tmpEpsilon);

            } else {

                if (missingCount == dimTrait) { // All missing, but not completely at random

                    simulateTraitForInternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);

                } else { // Partially observed

                    System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait); // copy observed values

                    final PartiallyMissingInformation.HashedIntArray indices = missingInformation.getMissingIndices(nodeIndex, traitIndex);
                    final int[] observed = indices.getComplement();
                    final int[] missing = indices.getArray();

                    final DenseMatrix64F V1 = new DenseMatrix64F(dimTrait, dimTrait);
                    CommonOps.scale(1.0 / branchPrecision, Vd, V1);

                    ConditionalVarianceAndTransform2 transform =
                            new ConditionalVarianceAndTransform2(
                                    V1, missing, observed
                            ); // TODO Cache (via delegated function)

                    final DenseMatrix64F cP0 = new DenseMatrix64F(missing.length, missing.length);
                    gatherRowsAndColumns(P0, cP0, missing, missing);

                    final WrappedVector cM2 = transform.getConditionalMean(
                            partialNodeBuffer, offsetPartial, // Tip value
                            sample, offsetParent); // Parent value

                    final DenseMatrix64F cP1 = transform.getConditionalPrecision();

                    final DenseMatrix64F cP2 = new DenseMatrix64F(missing.length, missing.length);
                    final DenseMatrix64F cV2 = new DenseMatrix64F(missing.length, missing.length);
                    CommonOps.add(cP0, cP1, cP2);

                    safeInvert(cP2, cV2, false);
                    double[][] cC2 = getCholeskyOfVariance(cV2.getData(), missing.length);

                    // TODO Drift?
                    assert (!likelihoodDelegate.getDiffusionProcessDelegate().hasDrift());

                    MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                            cM2, // input mean
                            cC2, 1.0, // input variance
                            new WrappedVector.Indexed(sample, offsetSample, missing, missing.length), // output sample
                            tmpEpsilon);


                    if (DEBUG) {
                        final WrappedVector M0 = new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait);

                        final WrappedVector M1 = new WrappedVector.Raw(sample, offsetParent, dimTrait);
                        final DenseMatrix64F P1 = new DenseMatrix64F(dimTrait, dimTrait);
                        CommonOps.scale(branchPrecision, Pd, P1);

                        final WrappedVector newSample = new WrappedVector.Raw(sample, offsetSample, dimTrait);

                        System.err.println("sT F E N");
                        System.err.println("M0: " + M0);
                        System.err.println("P0: " + P0);
                        System.err.println("");
                        System.err.println("M1: " + M1);
                        System.err.println("P1: " + P1);
                        System.err.println("");
                        System.err.println("cP0: " + cP0);
                        System.err.println("cM2: " + cM2);
                        System.err.println("cP1: " + cP1);
                        System.err.println("cP2: " + cP2);
                        System.err.println("cV2: " + cV2);
                        System.err.println("cC2: " + new Matrix(cC2));
                        System.err.println("SS: " + newSample);
                        System.err.println("");
                    }
                }
            }
        }
    }

    private WrappedVector getMeanWithDrift(double[] mean, int offsetMean, double[] drift, int dim) {
        for (int i = 0;i < dim; ++i) {
            tmpDrift[i] = mean[offsetMean + i] + drift[i];
        }
        return new WrappedVector.Raw(tmpDrift, 0, dimTrait);
    }

    private void simulateTraitForInternalNode(final int offsetSample,
                                              final int offsetParent,
                                              final int offsetPartial,
                                              final double branchPrecision) {

        if (!Double.isInfinite(branchPrecision)) {

            final WrappedVector M0 = new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait);
            final DenseMatrix64F P0 = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);

            final WrappedVector M1;
            final DenseMatrix64F P1;

            if (hasNoDrift) {
                M1 = new WrappedVector.Raw(sample, offsetParent, dimTrait);
                P1 = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.scale(branchPrecision, Pd, P1);
            } else {
                M1 = getMeanWithDrift(sample, offsetParent, displacementBuffer, dimTrait);
                P1 = DenseMatrix64F.wrap(dimTrait, dimTrait, precisionBuffer);
            }

//            boolean DEBUG_PRECISION = false;
//
//            if (DEBUG_PRECISION) {
//                DenseMatrix64F tP1 = new DenseMatrix64F(dimTrait, dimTrait);
//                CommonOps.scale(branchPrecision, Pd, tP1);
//
//                for (int i = 0; i < dimTrait; ++i) {
//                    for (int j = 0; j < dimTrait; ++j) {
//                        if (Math.abs(tP1.get(i,j) - P1.get(i,j)) != 0.0) {
//                            System.err.println("Unequal");
//                            System.exit(-1);
//                        }
//                    }
//                }
//            }

            final WrappedVector M2 = new WrappedVector.Raw(tmpMean, 0, dimTrait);
            final DenseMatrix64F P2 = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F V2 = new DenseMatrix64F(dimTrait, dimTrait);

            CommonOps.add(P0, P1, P2);
            safeInvert(P2, V2, false);
            weightedAverage(M0, P0, M1, P1, M2, V2, dimTrait);

            double[][] C2 = getCholeskyOfVariance(V2.getData(), dimTrait);

            MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                    M2.getBuffer(), 0, // input mean
                    C2, 1.0, // input variance
                    sample, offsetSample, // output sample
                    tmpEpsilon);

            if (DEBUG) {
                System.err.println("sT F I N");
                System.err.println("M0: " + M0);
                System.err.println("P0: " + P0);
                System.err.println("M1: " + M1);
                System.err.println("P1: " + P1);
                System.err.println("M2: " + M2);
                System.err.println("V2: " + V2);
                System.err.println("C2: " + new Matrix(C2));
                System.err.println("SS: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));
                System.err.println("");

                if (!check(M2))  {
                    System.exit(-1);
                }
            }

        } else {

            System.arraycopy(sample, offsetParent, sample, offsetSample, dimTrait);

            if (DEBUG) {
                System.err.println("sT F I N infinite branch precision");
                System.err.println("SS: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));
            }
        }
    }

    private boolean check(WrappedVector m2) {
        for (int i = 0; i < m2.getDim(); ++i) {
            if (Double.isNaN(m2.get(i))) {
                return false;
            }
        }
        return true;
    }
}
