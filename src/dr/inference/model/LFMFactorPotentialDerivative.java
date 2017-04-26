package dr.inference.model;

/**
 * Created by maxryandolinskytolkoff on 3/1/17.
 */
public class LFMFactorPotentialDerivative implements GradientProvider {
    LatentFactorModel lfm;

    public LFMFactorPotentialDerivative(LatentFactorModel lfm){
        this.lfm = lfm;
    }

    @Override
    public Likelihood getLikelihood() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getDimension() {
        return lfm.getFactors().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] derivative = new double[lfm.getFactors().getDimension()];
        Parameter missingIndicator = lfm.getMissingIndicator();
        int ntaxa = lfm.getFactors().getColumnDimension();
        int ntraits = lfm.getLoadings().getRowDimension();
        int nfac = lfm.getLoadings().getColumnDimension();
        double[] residual = lfm.getResidual();

        for (int i = 0; i < nfac; i++) {
            for (int j = 0; j < ntraits; j++) {
                for (int k = 0; k < ntaxa; k++) {
                    if(missingIndicator == null || missingIndicator.getParameterValue(k * ntraits + j) != 1){
                        derivative[k * nfac + i] -= lfm.getLoadings().getParameterValue(j, i) * lfm.getColumnPrecision().getParameterValue(j, j) *
                                residual[k * ntraits + j];
                    }
                }
            }
        }

        return derivative;
    }
}
