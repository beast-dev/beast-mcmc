package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

public interface SpeciationModelGradientProvider {

    default double getNodeHeightGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getBirthRateGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getDeathRateGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getSamplingRateGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getTreatmentProbabilityGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getSamplingProbabilityGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getBirthRateParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getDeathRateParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getSamplingRateParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getTreatmentProbabilityParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getSamplingProbabilityParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getBreakPoints() { throw new RuntimeException("Not yet implemented"); }

    // TODO This factory is probability unnecessary ...
    static SpeciationModelGradientProvider factory(SpeciationModel speciationModel) {
        return speciationModel.getProvider();
    }
}
