package dr.evomodel.treedatalikelihood.preorder;

import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

/**
 * @author Marc A. Suchard
 */
public class ConditionalVarianceAndTransform {

    /**
     * For partially observed tips: (y_1, y_2)^t \sim N(\mu, \Sigma) where
     * <p>
     * \mu = (\mu_1, \mu_2)^t
     * \Sigma = ((\Sigma_{11}, \Sigma_{12}), (\Sigma_{21}, \Sigma_{22})^t
     * <p>
     * then  y_1 | y_2 \sim N (\bar{\mu}, \bar{\Sigma}), where
     * <p>
     * \bar{\mu} = \mu_1 + \Sigma_{12}\Sigma_{22}^{-1}(y_2 - \mu_2), and
     * \bar{\Sigma} = \Sigma_{11} - \Sigma_{12}\Sigma_{22}^1\Sigma{21}
     */

    final private double[][] cholesky;
    final private Matrix affineTransform;
    private Matrix sBar;
    private final int[] missingIndices;
    private final int[] notMissingIndices;
    private final double[] tempStorage;

    private final int numMissing;
    private final int numNotMissing;

    private static final boolean DEBUG = false;

    public ConditionalVarianceAndTransform(final Matrix variance, final int[] missingIndices, final int[] notMissingIndices) {

        assert (missingIndices.length + notMissingIndices.length == variance.rows());
        assert (missingIndices.length + notMissingIndices.length == variance.columns());

        this.missingIndices = missingIndices;
        this.notMissingIndices = notMissingIndices;

        if (DEBUG) {
            System.err.println("variance:\n" + variance);
        }

        Matrix S12S22Inv = null;
        sBar = null;

        try {

            Matrix S22 = variance.extractRowsAndColumns(notMissingIndices, notMissingIndices);
            if (DEBUG) {
                System.err.println("S22:\n" + S22);
            }

            Matrix S22Inv = S22.inverse();
            if (DEBUG) {
                System.err.println("S22Inv:\n" + S22Inv);
            }

            Matrix S12 = variance.extractRowsAndColumns(missingIndices, notMissingIndices);
            if (DEBUG) {
                System.err.println("S12:\n" + S12);
            }

            S12S22Inv = S12.product(S22Inv);
            if (DEBUG) {
                System.err.println("S12S22Inv:\n" + S12S22Inv);
            }

            Matrix S12S22InvS21 = S12S22Inv.productWithTransposed(S12);
            if (DEBUG) {
                System.err.println("S12S22InvS21:\n" + S12S22InvS21);
            }

            sBar = variance.extractRowsAndColumns(missingIndices, missingIndices);
            sBar.decumulate(S12S22InvS21);
            if (DEBUG) {
                System.err.println("sBar:\n" + sBar);
            }

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        this.affineTransform = S12S22Inv;
        this.cholesky = ProcessSimulationDelegate.AbstractContinuousTraitDelegate.getCholeskyOfVariance(sBar);
        this.tempStorage = new double[missingIndices.length];

        this.numMissing = missingIndices.length;
        this.numNotMissing = notMissingIndices.length;

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
                delta += affineTransform.component(i, k) * shift[k];
            }

            muBar[i] = mu[offsetMu + missingIndices[i]] + delta;
        }

        return muBar;
    }

//    void scatterResult(final double[] source, final int offsetSource,
//                       final double[] destination, final int offsetDestination) {
//        for (int i = 0; i < numMissing; ++i) {
//            destination[offsetDestination + missingIndices[i]] = source[offsetSource + i];
//        }
//    }

    double[][] getConditionalCholesky() {
        return cholesky;
    }

    Matrix getVariance() {
        return sBar;
    }

//    Matrix getAffineTransform() {
//        return affineTransform;
//    }

    double[] getTemporaryStorage() {
        return tempStorage;
    }
}
