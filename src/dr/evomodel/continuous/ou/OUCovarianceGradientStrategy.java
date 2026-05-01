package dr.evomodel.continuous.ou;

public interface OUCovarianceGradientStrategy {

    void accumulate(double dt,
                    int dimension,
                    double[][] selection,
                    double[][] diffusion,
                    double[][] covarianceAdjoint,
                    double[] gradientAccumulator);

    void accumulateFlat(double dt,
                        int dimension,
                        double[][] selection,
                        double[][] diffusion,
                        double[] covarianceAdjoint,
                        boolean transposeAdjoint,
                        double[] gradientAccumulator);

    default void accumulateFlatMatrices(final double dt,
                                        final int dimension,
                                        final double[] selection,
                                        final double[] diffusion,
                                        final double[] covarianceAdjoint,
                                        final boolean transposeAdjoint,
                                        final double[] gradientAccumulator) {
        final double[][] denseSelection = new double[dimension][dimension];
        final double[][] denseDiffusion = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(selection, i * dimension, denseSelection[i], 0, dimension);
            System.arraycopy(diffusion, i * dimension, denseDiffusion[i], 0, dimension);
        }
        accumulateFlat(
                dt,
                dimension,
                denseSelection,
                denseDiffusion,
                covarianceAdjoint,
                transposeAdjoint,
                gradientAccumulator);
    }
}
