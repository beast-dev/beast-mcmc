package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

public class NewBirthDeathSerialSamplingModelGradient implements SpeciationModelGradientProvider {

    private final NewBirthDeathSerialSamplingModel model;

    public NewBirthDeathSerialSamplingModelGradient(NewBirthDeathSerialSamplingModel model) {
        this.model = model;
    }

    // Put gradient code here!
}
