package dr.evomodel.speciation;

// Incorporated into the base model

/**
 * The reconstructed birth-death model with incomplete sampling.
 */
//public class BirthDeathSamplingModel extends BirthDeathGernhard08Model {
//
//    public static final String BIRTH_DEATH_SAMPLING_MODEL = "birthDeathSamplingModel";
//
//    private Parameter samplingProportion;
//
//    public BirthDeathSamplingModel(
//            Parameter birthDiffRateParameter,
//            Parameter relativeDeathRateParameter,
//            Parameter samplingProportion,
//            Type units) {
//        this(BIRTH_DEATH_SAMPLING_MODEL, birthDiffRateParameter,
//                relativeDeathRateParameter, samplingProportion, units);
//    }
//
//
//    BirthDeathSamplingModel(String name,
//                            Parameter birthDiffRateParameter,
//                            Parameter relativeDeathRateParameter,
//                            Parameter samplingProportion,
//                            Type units) {
//
//        super(name, birthDiffRateParameter, relativeDeathRateParameter, units);
//
//        this.samplingProportion = samplingProportion;
//        addParameter(samplingProportion);
//        samplingProportion.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
//    }
//
//    public final double getS() {
//        return samplingProportion.getParameterValue(0);
//    }
//
//    public final double getLambda() {
//
//        double lambda = getR() / (1 - getA());
//
//        return lambda;
//    }
//
//    public final double getMu() {
//
//        double mu = getR() * getA() / (1 - getA());
//
//        return mu;
//    }
//
//
//    public double logTreeProbability(int taxonCount) {
//
//        int n = taxonCount;
//
//        return logGamma(n + 1) + (2 * n - 1) * Math.log(getR()) + Math.log(getLambda() * getS()) * (n - 1);
//
//    }
//
//    public double logNodeProbability(Tree tree, NodeRef node) {
//        final double height = tree.getNodeHeight(node);
//        final double mrh = -getR() * height;
//        final double z = Math.log(getS() * getLambda() + (getLambda() * (1 - getS()) - getMu()) * Math.exp(mrh));
//        double l = -2 * z + mrh;
//
//        if (tree.getRoot() == node) {
//            l += mrh - z;
//        }
//        return l;
//    }
//
//    public boolean includeExternalNodesInLikelihoodCalculation() {
//        return false;
//    }
//}
