package dr.evomodel.treelikelihood.thorneytreelikelihood;


public class BranchLengthLikelihoodDelegate implements ThorneyBranchLengthLikelihoodDelegate {
    private final double scale;

    public BranchLengthLikelihoodDelegate(double scale){
        this.scale = scale;
    }
    @Override
    public double getLogLikelihood(double observed, double expected) {
        return SaddlePointExpansion.logPoissonProbability(expected*scale, (int) Math.round(observed));
    }
}
