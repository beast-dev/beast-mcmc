package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

public class NewBirthDeathSerialSamplingModelGradient implements SpeciationModelGradientProvider {

    private final NewBirthDeathSerialSamplingModel model;

    public NewBirthDeathSerialSamplingModelGradient(NewBirthDeathSerialSamplingModel model) {
        this.model = model;
    }

    // Put gradient code here!
    private double partialQpartialRho(double t) {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();

        // d/dc2
        double partialC2 = 2 * lambda / constants[0];

        double partialQ = -4.0 * constants[1] * partialC2 + Math.exp(constants[0] * t )
                * -2.0 * (1 - constants[1]) * partialC2
                + Math.exp(constants[0] * t) * 2.0 * (1 + constants[1]) * partialC2;

        return partialQ;
    }

    public double[] getSamplingProbabilityGradient(Tree tree, NodeRef node) {
        model.precomputeConstants();

        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();

        // d/dc2
        double partialC2 = 2 * lambda / constants[0];

        // d/dq
        // TODO make right
        // TODO propagate Q properly
        double partialQ = 0.0;
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            double t = tree.getNodeHeight(tree.getInternalNode(i));
            partialQ += partialQpartialRho(t);
        }
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            double t = tree.getNodeHeight(tree.getExternalNode(i));
            partialQ += partialQpartialRho(t);
        }

        // TODO make right
        double gradientRho = partialC2 + partialQ;
        
        return new double[]{gradientRho};
    }


}
