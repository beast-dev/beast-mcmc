package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.distributions.MultivariateNormalDistribution;

import java.util.HashMap;

/**
 * @author Marc A. Suchard
 */
public class ConditionalOnPartiallyMissingTipsRealizedDelegate extends ConditionalOnTipsRealizedDelegate {

    public ConditionalOnPartiallyMissingTipsRealizedDelegate(String name, MutableTreeModel tree,
                                                             MultivariateDiffusionModel diffusionModel,
                                                             ContinuousTraitDataModel dataModel,
                                                             ConjugateRootTraitPrior rootPrior,
                                                             ContinuousRateTransformation rateTransformation,
                                                             BranchRateModel rateModel,
                                                             ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        this(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate,
                new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate));
    }

    public ConditionalOnPartiallyMissingTipsRealizedDelegate(String name, MutableTreeModel tree,
                                                             MultivariateDiffusionModel diffusionModel,
                                                             ContinuousTraitDataModel dataModel,
                                                             ConjugateRootTraitPrior rootPrior,
                                                             ContinuousRateTransformation rateTransformation,
                                                             BranchRateModel rateModel,
                                                             ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                             PartiallyMissingInformation missingInformation) {

        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
        this.missingInformation = missingInformation;

        assert (dataModel.getPrecisionType() == PrecisionType.FULL);

        throw new RuntimeException("Should remove this class");
    }

    @Override
    protected boolean isLoggable() {
        return false;
    }

    final private PartiallyMissingInformation missingInformation;

    @Override
    protected void simulateNode(final BranchNodeOperation operation,
                                final double branchNormalization) {
        final int nodeIndex = operation.getNodeNumber();
        likelihoodDelegate.getPostOrderPartial(nodeIndex, partialNodeBuffer);

        int offsetPartial = 0;
        int offsetSample = dimNode * nodeIndex;
        int offsetParent = dimNode * operation.getParentNumber();

        final boolean isExternal = nodeIndex < tree.getExternalNodeCount();

        final double branchPrecision = 1.0 / (operation.getBranchLength() * branchNormalization);

        for (int trait = 0; trait < numTraits; ++trait) {

            final double nodePrecision = partialNodeBuffer[offsetPartial + dimTrait];  // TODO PrecisionType.FULL

            if (!isExternal) {

                simulateTraitForNode(nodeIndex, trait, offsetSample, offsetParent, offsetPartial, nodePrecision);

            } else { // Is external

                // Copy tip values into sample
                System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);

                if (missingInformation.isPartiallyMissing(nodeIndex, trait)) {

                    PartiallyMissingInformation.HashedIntArray missingIndices =
                            missingInformation.getMissingIndices(nodeIndex, trait);

                    final int numMissing = missingIndices.getLength();
                    final int numNotMissing = missingIndices.getComplementLength();

                    assert (numMissing + numNotMissing == dimTrait);

                    ConditionalVarianceAndTranform transform;
                    try {
                        transform = conditionalMap.get(missingIndices);
                    } catch (NullPointerException nep) {
//                            System.err.println("Make CVT");
                        transform =
//                                    new ConditionalVarianceAndTranform(diffusionVariance,
//                                    missingIndices.getArray(),
//                                    missingIndices.getComplement());
                                null;

                        if (conditionalMap == null) {
                            conditionalMap = new HashMap<PartiallyMissingInformation.HashedIntArray,
                                    ConditionalVarianceAndTranform>();
                        }
                        conditionalMap.put(missingIndices, transform);
                    }
                    // TODO Must clear cache

//                        ConditionalVarianceAndTranform transform =
//                                new ConditionalVarianceAndTranform(diffusionVariance,
//                                        missingIndices.getArray(),
//                                        missingIndices.getComplement());

                    // TODO PrecisionType.FULL

                    final double[] conditionalMean = transform.getConditionalMean(sample, offsetSample,
                            sample, offsetParent);
                    final double[][] conditionalCholesky = transform.getConditionalCholesky();

                    final double sqrtScale = Math.sqrt(1.0 / branchPrecision);

                    MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                            conditionalMean, 0, // input mean
                            conditionalCholesky, sqrtScale, // input variance
                            tmpMean, 0, // output sample
                            transform.getTemporageStorage());

                    for (int i = 0; i < numMissing; ++i) {
                        sample[offsetSample + missingIndices.get(i)] = tmpMean[i];
                    }

//                        System.err.println("mean:\n" + new Vector(conditionalMean));
//                        System.err.println("cholesky:\n" + new Matrix(conditionalCholesky));
//                        System.err.println("sS: " + sqrtScale);
//                        System.err.println("cMean\n" + new Vector(tmpMean));
//                        System.err.println("");
//                        System.err.println("");
                }
            }

            offsetSample += dimTrait;
            offsetParent += dimTrait;
            offsetPartial += (dimTrait + 1); // TODO PrecisionType.FULL
        }
    }

    public static final String PARTIAL = "partial";

    public static String getPartiallyMissingTraitName(final String traitName) {
        return PARTIAL + "." + traitName;
    }
}
