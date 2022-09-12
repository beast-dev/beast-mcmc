package dr.inference.operators.repeatedMeasures;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ConditionalTraitSimulationHelper;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Attribute;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.getTipTraitName;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */
public interface GammaGibbsProvider {


    SufficientStatistics getSufficientStatistics(int dim);

    Parameter getPrecisionParameter();

    void drawValues();

    class SufficientStatistics {
        final public int observationCount;
        final public double sumOfSquaredErrors;

        public SufficientStatistics(int observationCount, double sumOfSquaredErrors) {
            this.observationCount = observationCount;
            this.sumOfSquaredErrors = sumOfSquaredErrors;
        }
    }

    class Default implements GammaGibbsProvider {

        private final Parameter precisionParameter;
        private final Parameter meanParameter;
        private final boolean isLog;
        private final List<Attribute<double[]>> dataList;

        public Default(DistributionLikelihood inLikelihood) {

            Distribution likelihood = inLikelihood.getDistribution();
            this.dataList = inLikelihood.getDataList();

            if (likelihood instanceof NormalDistributionModel) {
                this.precisionParameter = (Parameter) ((NormalDistributionModel) likelihood).getPrecision();
                this.meanParameter = (Parameter) ((NormalDistributionModel) likelihood).getMean();
                this.isLog = false;
            } else if (likelihood instanceof LogNormalDistributionModel) {
                if (((LogNormalDistributionModel) likelihood).getParameterization() == LogNormalDistributionModel.Parameterization.MU_PRECISION) {
                    this.meanParameter = ((LogNormalDistributionModel) likelihood).getMuParameter();
                } else {
                    throw new RuntimeException("Must characterize likelihood in terms of mu and precision parameters");
                }
                this.precisionParameter = ((LogNormalDistributionModel) likelihood).getPrecisionParameter();
                isLog = true;
            } else
                throw new RuntimeException("Likelihood must be Normal or log Normal");

            if (precisionParameter == null)
                throw new RuntimeException("Must characterize likelihood in terms of a precision parameter");
        }

        @Override
        public SufficientStatistics getSufficientStatistics(int dim) {

            // Calculate weighted sum-of-squares
            final double mu = meanParameter.getParameterValue(dim);
            double SSE = 0;
            int n = 0;
            for (Attribute<double[]> statistic : dataList) {
                for (double x : statistic.getAttributeValue()) {
                    if (isLog) {
                        final double logX = Math.log(x);
                        SSE += (logX - mu) * (logX - mu);
                    } else {
                        SSE += (x - mu) * (x - mu);
                    }
                    n++;
                }
            }

            return new SufficientStatistics(n, SSE);
        }

        @Override
        public Parameter getPrecisionParameter() {
            return precisionParameter;
        }

        @Override
        public void drawValues() {
            // Do nothing
        }
    }

    class NormalExtensionGibbsProvider implements GammaGibbsProvider {

        private final ModelExtensionProvider.NormalExtensionProvider dataModel;
        private final TreeDataLikelihood treeLikelihood;
        private final ConditionalTraitSimulationHelper traitProvider;
        private final TreeTrait tipTrait;
        private final CompoundParameter traitParameter;
        private final Parameter precisionParameter;
        private final boolean[] missingVector;

        private double[] tipValues;
        private boolean hasCheckedDimension = false;

        public NormalExtensionGibbsProvider(ModelExtensionProvider.NormalExtensionProvider dataModel,
                                            TreeDataLikelihood treeLikelihood,
                                            String traitName) {
            this.dataModel = dataModel;
            this.treeLikelihood = treeLikelihood;
            this.traitParameter = dataModel.getParameter();
            this.missingVector = dataModel.getDataMissingIndicators();

            if (traitName == null) {
                this.tipTrait = null;
                this.traitProvider = ((ContinuousDataLikelihoodDelegate)
                        treeLikelihood.getDataLikelihoodDelegate()).getExtensionHelper();
            } else {
                this.tipTrait = treeLikelihood.getTreeTrait(getTipTraitName(traitName));
                this.traitProvider = null;
            }


            MatrixParameterInterface matrixParameter = dataModel.getExtensionPrecision();

            if (matrixParameter instanceof DiagonalMatrix) {
                this.precisionParameter = ((DiagonalMatrix) matrixParameter).getDiagonalParameter();

            } else { //TODO: alternatively, check that the off-diagonal elements are zero every time you update the parameter?
                //TODO: does this belong in the parser?
                throw new RuntimeException(this.getClass().getName() +
                        " only applies to diagonal precision matrices, but the " +
                        ModelExtensionProvider.NormalExtensionProvider.class.getName() +
                        " supplied a precision matrix of class " +
                        matrixParameter.getClass().getName() + ".");
            }

        }

        @Override
        public SufficientStatistics getSufficientStatistics(int dim) {

            final int taxonCount = treeLikelihood.getTree().getExternalNodeCount();
            final int traitDim = dataModel.getDataDimension();

            if (!hasCheckedDimension) { //TODO: actually check that this works
                if (taxonCount * traitDim != tipValues.length) {
                    throw new RuntimeException("dimensions are incompatible");
                }
                hasCheckedDimension = true;
            }

            int missingCount = 0;

            double SSE = 0;

            for (int taxon = 0; taxon < taxonCount; ++taxon) {

                int offset = traitDim * taxon;
                if (missingVector == null || !missingVector[dim + offset]) {
                    double traitValue = traitParameter.getParameter(taxon).getParameterValue(dim);
                    double tipValue = tipValues[taxon * traitDim + dim];

                    SSE += (traitValue - tipValue) * (traitValue - tipValue);
                } else {
                    missingCount += 1;
                }
            }

            return new SufficientStatistics(taxonCount - missingCount, SSE);
        }

        @Override
        public Parameter getPrecisionParameter() {
            return precisionParameter;
        }

        @Override
        public void drawValues() {
            double[] tipTraits = tipTrait == null ? traitProvider.drawTraitsAbove(dataModel) :
                    (double[]) tipTrait.getTrait(treeLikelihood.getTree(), null);
            tipValues = dataModel.transformTreeTraits(tipTraits);
            if (DEBUG) {
                System.err.println("tipValues: " + new WrappedVector.Raw(tipValues));
            }
        }

        private static final boolean DEBUG = false;
    }

    class MultiplicativeGammaGibbsProvider implements GammaGibbsProvider {
        // TODO: add citation to "Sparse Bayesian infinite factor models BY A. BHATTACHARYA AND D. B. DUNSON, Biometrika (2011)"

        private final Parameter rowMultipliers;
        private final MultiplicativeGammaGibbsHelper helper;

        public MultiplicativeGammaGibbsProvider(Parameter rowMultipliers,
                                                MultiplicativeGammaGibbsHelper helper) {
            this.rowMultipliers = rowMultipliers;
            this.helper = helper;
        }


        @Override
        public SufficientStatistics getSufficientStatistics(int dim) {


            double rateSum = 0;

            int k = helper.getColumnDimension();
            int p = helper.getRowDimension();

            for (int i = dim; i < k; i++) {
                double globalConst = gpMult(i + 1, dim);
                double sumSquares = helper.computeSumSquaredErrors(i);

                rateSum += globalConst * sumSquares;
            }

            return new SufficientStatistics(p * (k - dim), rateSum);

        }

        @Override
        public Parameter getPrecisionParameter() {
            return rowMultipliers;
        }

        @Override
        public void drawValues() {
            // do nothing
        }

        private double gpMult(int multTo, int skip) { // TODO: probably could be more efficient
            double value = 1.0;
            for (int i = 0; i < multTo; i++) {
                if (i != skip) { // TODO: could remove 'if' statement with two for loops (probably doesn't matter)
                    value *= rowMultipliers.getParameterValue(i);
                }
            }
            return value;
        }

    }


}
