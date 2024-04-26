package dr.evomodel.treedatalikelihood.preorder;

import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

public class ConditionalPrecisionAndTransform2 {

    /**
     * For partially observed tips: (y_1, y_2)^t \sim N(\mu, P^{-1}) where
     * <p>
     * \mu = (\mu_1, \mu_2)^t
     * \Sigma = ((P_{11}, P_{12}), (P_{21}, P_{22})^t
     * <p>
     * then  y_1 | y_2 \sim N (\bar{\mu}, \bar{P}^{-1}), where
     * <p>
     * \bar{\mu} = \mu_1 - P_{11}^{-1}P_{12}(y_2 - \mu_2), and
     * \bar{P} = P_{11}
     */

    private final DenseMatrix64F affineTransform;
    private final DenseMatrix64F P11Inv;

    private final int[] missingIndices;
    private final int[] notMissingIndices;

    private final int numMissing;
    private final int numNotMissing;

    private static final boolean DEBUG = false;

    public ConditionalPrecisionAndTransform2(final DenseMatrix64F precision,
                                             final int[] missingIndices,
                                             final int[] notMissingIndices) {

        assert (missingIndices.length + notMissingIndices.length == precision.getNumRows());
        assert (missingIndices.length + notMissingIndices.length == precision.getNumCols());

        this.missingIndices = missingIndices;
        this.notMissingIndices = notMissingIndices;

        DenseMatrix64F P11 = new DenseMatrix64F(missingIndices.length, missingIndices.length);
        MissingOps.gatherRowsAndColumns(precision, P11, missingIndices, missingIndices);

        P11Inv = new DenseMatrix64F(missingIndices.length, missingIndices.length);
        CommonOps.invert(P11, P11Inv);

        DenseMatrix64F P12 = new DenseMatrix64F(missingIndices.length, notMissingIndices.length);
        MissingOps.gatherRowsAndColumns(precision, P12, missingIndices, notMissingIndices);

        DenseMatrix64F P11InvP12 = new DenseMatrix64F(missingIndices.length, notMissingIndices.length);
        CommonOps.mult(P11Inv, P12, P11InvP12);

        this.affineTransform = P11InvP12;

        this.numMissing = missingIndices.length;
        this.numNotMissing = notMissingIndices.length;
    }

    public DenseMatrix64F getConditionalVariance() {
        return P11Inv;
    }

    public double[] getConditionalMean(final double[] y, final int offsetY,
                                       final double[] mu, final int offsetMu) {

        double[] muBar = new double[numMissing];

        double[] shift = new double[numNotMissing];
        for (int i = 0; i < numNotMissing; ++i) {
            final int notI = notMissingIndices[i];
            shift[i] = y[offsetY + notI] - mu[offsetMu + notI];
        }

        for (int i = 0; i < numMissing; ++i) {
            double delta = 0.0;
            for (int k = 0; k < numNotMissing; ++k) {
                delta += affineTransform.unsafe_get(i, k) * shift[k];
            }

            muBar[i] = mu[offsetMu + missingIndices[i]] - delta;
        }

        return muBar;
    }
}

