package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import org.ejml.data.DenseMatrix64F;

final class OrthogonalBlockCanonicalTransitionAssembler {

    private OrthogonalBlockCanonicalTransitionAssembler() { }

    static void fillCanonicalTransition(final DenseMatrix64F transitionMatrix,
                                        final DenseMatrix64F transitionCovariance,
                                        final double[] stationaryMean,
                                        final double[] transitionOffsetScratch,
                                        final double[] transitionCovarianceScratch,
                                        final double[] choleskyScratch,
                                        final double[] lowerInverseScratch,
                                        final CanonicalGaussianTransition out) {
        final int dimension = transitionMatrix.numRows;
        final double[] transitionData = transitionMatrix.data;
        fillTransitionOffset(transitionMatrix, stationaryMean, transitionOffsetScratch);

        final double logDet = OrthogonalBlockPositiveDefiniteInverter.copyAndInvertFlat(
                transitionCovariance,
                transitionCovarianceScratch,
                out.precisionYY,
                choleskyScratch,
                lowerInverseScratch);

        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += out.precisionYY[iOff + k] * transitionData[k * dimension + j];
                }
                out.precisionYX[iOff + j] = -sum;
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out.precisionXY[i * dimension + j] = out.precisionYX[j * dimension + i];
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum -= transitionData[k * dimension + i] * out.precisionYX[k * dimension + j];
                }
                out.precisionXX[i * dimension + j] = sum;
            }
        }

        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                sum += out.precisionYY[iOff + j] * transitionOffsetScratch[j];
            }
            out.informationY[i] = sum;
        }

        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimension; ++j) {
                sum += transitionData[j * dimension + i] * out.informationY[j];
            }
            out.informationX[i] = -sum;
        }

        double quadratic = 0.0;
        for (int i = 0; i < dimension; ++i) {
            quadratic += transitionOffsetScratch[i] * out.informationY[i];
        }
        out.logNormalizer = 0.5 * (dimension * Math.log(2.0 * Math.PI) + logDet + quadratic);
    }

    static void fillTransitionOffset(final DenseMatrix64F transitionMatrix,
                                     final double[] stationaryMean,
                                     final double[] out) {
        final int dimension = transitionMatrix.numRows;
        final double[] transitionData = transitionMatrix.data;
        for (int i = 0; i < dimension; ++i) {
            double transformedMean = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                transformedMean += transitionData[rowOffset + j] * stationaryMean[j];
            }
            out[i] = stationaryMean[i] - transformedMean;
        }
    }
}
