package dr.inference.operators.repeatedMeasures;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Attribute;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

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

        SufficientStatistics(int observationCount, double sumOfSquaredErrors) {
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

    class RepeatedMeasuresGibbsProvider implements GammaGibbsProvider {

        //        private final RepeatedMeasuresTraitDataModel dataModel;
        private final TreeDataLikelihood treeLikelihood;
        private final CompoundParameter traitParameter;
        private final Parameter precisionParameter;
        private final TreeTrait tipTrait;
        private final boolean[] missingVector;

        private double tipValues[];

        public RepeatedMeasuresGibbsProvider(RepeatedMeasuresTraitDataModel dataModel,
                                             TreeDataLikelihood treeLikelihood,
                                             String traitName) {
//            this.dataModel = dataModel;
            this.treeLikelihood = treeLikelihood;
            this.traitParameter = dataModel.getParameter();
            this.precisionParameter = dataModel.getSamplingPrecision();
            this.tipTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);
            this.missingVector = dataModel.getMissingIndicator();
        }

        @Override
        public SufficientStatistics getSufficientStatistics(int dim) {

            final int taxonCount = treeLikelihood.getTree().getExternalNodeCount();
            final int traitDim = treeLikelihood.getDataLikelihoodDelegate().getTraitDim();
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
            tipValues = (double[]) tipTrait.getTrait(treeLikelihood.getTree(), null);
            if (DEBUG) {
                System.err.println("tipValues: " + new WrappedVector.Raw(tipValues));
            }
        }

        private static final boolean DEBUG = false;
    }

    //    TODO: Eliminate code duplication by merging FactorGibbsProvider with RepeatedMeasuresGibbsProvider
    class FactorGibbsProvider implements GammaGibbsProvider {

        private static final Boolean DEBUG = false;


        private final TreeDataLikelihood treeLikelihood;
        private final CompoundParameter traitParameter;
        private final Parameter precisionParameter;
        private final TreeTrait tipTrait;
        private final boolean[] missingVector;

        private final double tipValues[];

        private final IntegratedFactorAnalysisLikelihood dataModel;


        public FactorGibbsProvider(IntegratedFactorAnalysisLikelihood dataModel, TreeDataLikelihood treeLikelihood,
                                   String traitName) {

            this.dataModel = dataModel;
            this.treeLikelihood = treeLikelihood;
            this.traitParameter = dataModel.getParameter();
            this.precisionParameter = dataModel.getPrecision();
            this.tipTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);
            this.missingVector = dataModel.getMissingIndicator();

            int dimTaxa = dataModel.getNumberOfTaxa();
            int dimData = dataModel.getDataDimension();

            this.tipValues = new double[dimData * dimTaxa];


        }


        @Override
        public SufficientStatistics getSufficientStatistics(int dim) {
            final int taxonCount = treeLikelihood.getTree().getExternalNodeCount();
            final int traitDim = dataModel.getDataDimension();

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

            final int treeTraitDim = treeLikelihood.getDataLikelihoodDelegate().getTraitDim();
            final int traitDim = dataModel.getDataDimension();
            final int taxonCount = treeLikelihood.getTree().getExternalNodeCount();

            double[] treeTipValues = (double[]) tipTrait.getTrait(treeLikelihood.getTree(), null);
            MatrixParameterInterface loadings = dataModel.getLoadings();

            DenseMatrix64F treeTraitMatrix = DenseMatrix64F.wrap(taxonCount, treeTraitDim, treeTipValues);
            DenseMatrix64F loadingsMatrix = DenseMatrix64F.wrap(treeTraitDim, traitDim, loadings.getParameterValues());


            DenseMatrix64F traitMatrix = DenseMatrix64F.wrap(taxonCount, traitDim, tipValues);
            org.ejml.ops.CommonOps.mult(treeTraitMatrix, loadingsMatrix, traitMatrix);

            if (DEBUG) {
                treeTraitMatrix.print();
                loadingsMatrix.print();
                traitMatrix.print();
            }

        }
    }
}
