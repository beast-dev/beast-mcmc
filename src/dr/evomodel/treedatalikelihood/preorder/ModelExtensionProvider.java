package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;

public interface ModelExtensionProvider extends ContinuousTraitPartialsProvider{

    AbstractContinuousExtensionDelegate getExtensionDelegate(
            ProcessSimulationDelegate.AbstractContinuousTraitDelegate processDelegate,
            String traitName);

}
