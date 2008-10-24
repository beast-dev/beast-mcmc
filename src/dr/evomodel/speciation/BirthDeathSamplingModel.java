package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * The reconstructed birth-death model with incomplete sampling.
 */
public class BirthDeathSamplingModel extends BirthDeathGernhard08Model {

    public static final String BIRTH_DEATH_SAMPLING_MODEL = "birthDeathSamplingModel";

    private Parameter samplingProportion;

    public BirthDeathSamplingModel(
            Parameter birthDiffRateParameter,
            Parameter relativeDeathRateParameter,
            Parameter samplingProportion,
            Type units) {

        super(BIRTH_DEATH_SAMPLING_MODEL, birthDiffRateParameter, relativeDeathRateParameter, units);

        this.samplingProportion = samplingProportion;
        addParameter(samplingProportion);
        samplingProportion.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
    }

    public final double getS() {
        return samplingProportion.getParameterValue(0);
    }

    public double logTreeProbability(int taxonCount) {
        throw new RuntimeException("Not implemented yet.");
    }

    public double logNodeProbability(Tree tree, NodeRef node) {
        throw new RuntimeException("Not implemented yet.");
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return false;
    }
}
