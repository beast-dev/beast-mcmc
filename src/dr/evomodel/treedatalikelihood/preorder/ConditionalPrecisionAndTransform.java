package dr.evomodel.treedatalikelihood.preorder;

import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;


public class ConditionalPrecisionAndTransform {

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
    //final private double[][] cholesky;
    final private Matrix affineTransform;
    private Matrix pBar;
    private final int[] missingIndices;
    private final int[] notMissingIndices;

    private final int numMissing;
    private final int numNotMissing;

    private static final boolean DEBUG = false;

    public ConditionalPrecisionAndTransform(final Matrix precision, final int[] missingIndices, final int[] notMissingIndices) {

        assert (missingIndices.length + notMissingIndices.length == precision.rows());
        assert (missingIndices.length + notMissingIndices.length == precision.columns());

        this.missingIndices = missingIndices;
        this.notMissingIndices = notMissingIndices;

        if (DEBUG) {
            System.err.println("variance:\n" + precision);
        }

        Matrix P11InvP12 = null;
        pBar = null;

        try {
            Matrix P11 = precision.extractRowsAndColumns(missingIndices, missingIndices);
            if (DEBUG) {
                System.err.println("P11:\n" + P11);
            }

            Matrix P11Inv = P11.inverse();
            if (DEBUG) {
                System.err.println("P11Inv:\n" + P11Inv);
            }

            Matrix P12 = precision.extractRowsAndColumns(missingIndices, notMissingIndices);
            if (DEBUG) {
                System.err.println("P12:\n" + P12);
            }

            P11InvP12 = P11Inv.product(P12);
            if (DEBUG) {
                System.err.println("P11InvP12:\n" + P11InvP12);
            }

            pBar = precision.extractRowsAndColumns(missingIndices, missingIndices);
            if (DEBUG) {
                System.err.println("pBar:\n" + pBar);
            }


        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        this.affineTransform = P11InvP12;
        //this.cholesky = ProcessSimulationDelegate.AbstractContinuousTraitDelegate.getCholeskyOfVariance(sBar);
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

            muBar[i] = mu[offsetMu + missingIndices[i]] - delta;
        }

        return muBar;
    }


    public Matrix getConditionalPrecision() {
        return pBar;
    }

}

