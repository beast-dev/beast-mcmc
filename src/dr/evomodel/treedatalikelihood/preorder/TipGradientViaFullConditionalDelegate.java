package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.MatrixParameterInterface;

/**
 * @author Marc A. Suchard
 */
public class TipGradientViaFullConditionalDelegate extends TipFullConditionalDistributionDelegate {

    public TipGradientViaFullConditionalDelegate(String name, MutableTreeModel tree,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 ContinuousTraitDataModel dataModel,
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

        if (likelihoodDelegate.getPrecisionType() != PrecisionType.SCALAR) {
            throw new RuntimeException("Tip gradients are not implemented for '" +
                    likelihoodDelegate.getPrecisionType().toString() + "' likelihoods");
        }

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
                        precision.getParameterValue(i * dimTrait + j);
            }

            gradient[i] = sum;
        }

        return gradient;
    }
}
