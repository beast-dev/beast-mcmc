package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

/**
 * @author Marc A. Suchard
 */
public class TipGradientViaFullConditionalDelegate extends TipFullConditionalDistributionDelegate {

    public TipGradientViaFullConditionalDelegate(String name, MutableTreeModel tree,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 ContinuousTraitPartialsProvider dataModel,
                                                 ConjugateRootTraitPrior rootPrior,
                                                 ContinuousRateTransformation rateTransformation,
                                                 ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
    }

    public static String getName(String name) {
        return "grad." + name;
    }

    public String getTraitName(String name) {
        return getName(name);
    }

    @Override
    protected double[] getTraitForNode(NodeRef node) {

        if (likelihoodDelegate.getPrecisionType() == PrecisionType.SCALAR) {
            return getTraitForNodeScalar(node);
        } else if (likelihoodDelegate.getPrecisionType() == PrecisionType.FULL) {
            return getTraitForNodeFull(node);
        } else {
            throw new RuntimeException("Tip gradients are not implemented for '" +
                    likelihoodDelegate.getPrecisionType().toString() + "' likelihoods");
        }
    }

    private double[] getTraitForNodeScalar(NodeRef node) {

        final double[] fullConditionalPartial = super.getTraitForNode(node);

        final double[] postOrderPartial = new double[dimPartial * numTraits];
        cdi.getPostOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), postOrderPartial);

        final MatrixParameterInterface precision = diffusionModel.getPrecisionParameter();

        final double[] gradient = new double[dimTrait * numTraits];

        if (numTraits > 1) {
            throw new RuntimeException("Not yet implemented");
        }

        final double scale = fullConditionalPartial[dimTrait];

        for (int i = 0; i < dimTrait; ++i) {

            double sum = 0.0;
            for (int j = 0; j < dimTrait; ++j) {
                sum += (fullConditionalPartial[j] - postOrderPartial[j]) * scale *
                        precision.getParameterValue(i, j);
//                        precision.getParameterValue(i * dimTrait + j);
            }

            gradient[i] = sum;
        }

        return gradient;
    }

    protected double[] getTraitForNodeFull(NodeRef node) {

        if (numTraits > 1) {
            throw new RuntimeException("Not yet implemented");
        }

        // Exclude data models that have not been tested.
        if ((dataModel instanceof JointPartialsProvider)
                || (dataModel instanceof IntegratedProcessTraitDataModel)
                || (dataModel instanceof ElementaryVectorDataModel)
                || (dataModel instanceof EmptyTraitDataModel)) {
            throw new RuntimeException("Tip gradients are not implemented for '"
                    + dataModel.getClass().toString() + "' data model.");
        }

        // Pre stats
        final double[] fullConditionalPartial = super.getTraitForNode(node);
        NormalSufficientStatistics statPre = new NormalSufficientStatistics(fullConditionalPartial, 0, dimTrait, Pd, likelihoodDelegate.getPrecisionType());
        DenseMatrix64F precisionPre = statPre.getRawPrecisionCopy();

        if (dataModel instanceof ModelExtensionProvider.NormalExtensionProvider) {
            // Exclude data models that have not been tested.
            if (!(dataModel instanceof RepeatedMeasuresTraitDataModel)) {
                throw new RuntimeException("Tip gradients are not implemented for '"
                        + dataModel.getClass().toString() + "' data model.");
            }
            // TODO Move to different class ?
            // TODO Only works for one data point (no repetition) ?
            DenseMatrix64F samplingVariance = ((ModelExtensionProvider.NormalExtensionProvider) dataModel).getExtensionVariance(node);
            CommonOps.addEquals(samplingVariance, statPre.getRawVariance());
            safeInvert2(samplingVariance, precisionPre, false);
        }

        // Post mean
        final double[] postOrderPartial = new double[dimPartial * numTraits];
        int nodeIndex = likelihoodDelegate.getActiveNodeIndex(node.getNumber());
        cdi.getPostOrderPartial(nodeIndex, postOrderPartial);
        DenseMatrix64F meanPost = MissingOps.wrap(postOrderPartial, 0, dimTrait, 1);

        // - Q_i * (X_i - m_i)
        DenseMatrix64F gradient = new DenseMatrix64F(dimTrait, numTraits);
        CommonOps.addEquals(meanPost, -1.0, statPre.getRawMean());
        CommonOps.changeSign(meanPost);
        CommonOps.mult(precisionPre, meanPost, gradient);

        return gradient.getData();
    }
}
