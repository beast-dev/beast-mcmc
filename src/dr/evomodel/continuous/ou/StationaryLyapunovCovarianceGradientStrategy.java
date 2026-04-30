package dr.evomodel.continuous.ou;

final class StationaryLyapunovCovarianceGradientStrategy implements OUCovarianceGradientStrategy {

    @Override
    public void accumulate(final OUProcessModel model,
                           final double dt,
                           final int dimension,
                           final double[][] selection,
                           final double[][] diffusion,
                           final double[][] covarianceAdjoint,
                           final double[] gradientAccumulator) {
        model.accumulateViaStationaryLyapunov(
                dt, dimension, selection, diffusion, covarianceAdjoint, gradientAccumulator);
    }

    @Override
    public void accumulateFlat(final OUProcessModel model,
                               final double dt,
                               final int dimension,
                               final double[][] selection,
                               final double[][] diffusion,
                               final double[] covarianceAdjoint,
                               final boolean transposeAdjoint,
                               final double[] gradientAccumulator) {
        model.accumulateViaStationaryLyapunovFlat(
                dt, dimension, selection, diffusion, covarianceAdjoint, transposeAdjoint, gradientAccumulator);
    }
}
