package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;

public interface ModelExtensionProvider {

    ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate,
                                                     TreeTrait treeTrait,
                                                     Tree tree);

    interface NormalExtensionProvider extends ModelExtensionProvider, ContinuousTraitPartialsProvider {

        DenseMatrix64F getExtensionVariance();

        MatrixParameterInterface getExtensionPrecision();

        double[] transformTreeTraits(double[] treeTraits);
    }
}


