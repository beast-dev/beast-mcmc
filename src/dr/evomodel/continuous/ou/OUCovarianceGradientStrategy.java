package dr.evomodel.continuous.ou;

public interface OUCovarianceGradientStrategy {

    void accumulateFlat(double dt,
                        int dimension,
                        double[] selection,
                        double[] diffusion,
                        double[] covarianceAdjoint,
                        boolean transposeAdjoint,
                        double[] gradientAccumulator);
}
