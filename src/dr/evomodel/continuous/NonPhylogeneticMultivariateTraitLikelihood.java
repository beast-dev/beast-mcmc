package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import java.util.List;

/**
 * Integrated multivariate trait likelihood that assumes a fully-conjugate prior on the root and
 * no underlying tree structure.
 *
 * @author Gabriela Cybis
 * @author Marc A. Suchard
 */
public class NonPhylogeneticMultivariateTraitLikelihood extends FullyConjugateMultivariateTraitLikelihood {

    public NonPhylogeneticMultivariateTraitLikelihood(String traitName,
                                                     TreeModel treeModel,
                                                     MultivariateDiffusionModel diffusionModel,
                                                     CompoundParameter traitParameter,
                                                     Parameter deltaParameter,
                                                     List<Integer> missingIndices,
                                                     boolean cacheBranches,
                                                     boolean scaleByTime,
                                                     boolean useTreeLength,
                                                     BranchRateModel rateModel,
                                                     Model samplingDensity,
                                                     boolean reportAsMultivariate,
                                                     double[] rootPriorMean,
                                                     double rootPriorSampleSize,
                                                     boolean reciprocalRates) {
        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches,
                scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate, rootPriorMean,
                rootPriorSampleSize, reciprocalRates);
    }
   
    void postOrderTraverse(TreeModel treeModel, NodeRef node, double[][] precisionMatrix,
                           double logDetPrecisionMatrix, boolean cacheOuterProducts) {
        // TODO Need to write this function to avoid peeling on the tree
    }

    private static final boolean DEBUG_NO_TREE = true;
}