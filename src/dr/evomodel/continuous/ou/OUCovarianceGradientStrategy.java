package dr.evomodel.continuous.ou;

interface OUCovarianceGradientStrategy {

    void accumulate(OUProcessModel model,
                    double dt,
                    int dimension,
                    double[][] selection,
                    double[][] diffusion,
                    double[][] covarianceAdjoint,
                    double[] gradientAccumulator);

    void accumulateFlat(OUProcessModel model,
                        double dt,
                        int dimension,
                        double[][] selection,
                        double[][] diffusion,
                        double[] covarianceAdjoint,
                        boolean transposeAdjoint,
                        double[] gradientAccumulator);
}
