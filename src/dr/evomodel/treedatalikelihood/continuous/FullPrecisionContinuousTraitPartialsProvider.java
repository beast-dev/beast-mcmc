package dr.evomodel.treedatalikelihood.continuous;

import dr.inference.model.MatrixParameterInterface;

public interface FullPrecisionContinuousTraitPartialsProvider extends ContinuousTraitPartialsProvider {

    MatrixParameterInterface getExtensionPrecisionParameter();

}
