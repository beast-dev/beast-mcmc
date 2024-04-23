package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

public interface ModelExtensionProvider extends ContinuousTraitPartialsProvider {

    ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate,
                                                     TreeTrait treeTrait,
                                                     Tree tree);


    interface NormalExtensionProvider extends ModelExtensionProvider {

        boolean diagonalVariance();

        DenseMatrix64F getExtensionVariance();

        DenseMatrix64F getExtensionVariance(NodeRef node);

        MatrixParameterInterface getExtensionPrecision();

        void chainRuleWrtVariance(double[] gradient, NodeRef node);

        @Override
        default void updateTipDataGradient(DenseMatrix64F precision, DenseMatrix64F variance,
                                           NodeRef node, int offset, int dimGradient) {

            extendTipDataGradient(this, precision, variance, node, offset, dimGradient);
        }

        static void extendTipDataGradient(NormalExtensionProvider provider,
                                          DenseMatrix64F precision, DenseMatrix64F variance,
                                          NodeRef node, int offset, int dimGradient) {

            if (offset != 0 || dimGradient != provider.getTraitDimension()) {
                throw new RuntimeException("not implemented for subset of model.");
            }

            DenseMatrix64F samplingVariance = provider.getExtensionVariance(node);
            CommonOps.addEquals(samplingVariance, variance);
            safeInvert2(samplingVariance, precision, false);
        }

        @Override
        default boolean needToUpdateTipDataGradient(int offset, int dimGradient) {
            return true;
        }

    }
}


