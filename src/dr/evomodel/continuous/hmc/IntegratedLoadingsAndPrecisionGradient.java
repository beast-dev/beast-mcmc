package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.evomodel.treedatalikelihood.preorder.WrappedNormalSufficientStatistics;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.util.TaskPool;

public class IntegratedLoadingsAndPrecisionGradient extends IntegratedLoadingsGradient {

    CompoundParameter jointParameter;

    public IntegratedLoadingsAndPrecisionGradient(CompoundParameter jointParameter,
                                                  TreeDataLikelihood treeDataLikelihood,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                  IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                                  ContinuousTraitPartialsProvider partialsProvider,
                                                  TaskPool taskPool,
                                                  ThreadUseProvider threadUseProvider,
                                                  RemainderCompProvider remainderCompProvider) {
        super(treeDataLikelihood, likelihoodDelegate, factorAnalysisLikelihood, partialsProvider, taskPool, threadUseProvider, remainderCompProvider);
        this.jointParameter = jointParameter;
    }

    @Override
    public Parameter getParameter() {
        return jointParameter;
    }

    @Override
    protected int getGradientDimension() {
        return dimFactors * dimTrait + dimTrait;
    }

    @Override
    public int getDimension() {
        return getGradientDimension();
    }

    private void computePrecisionGradientForOneTaxon(int index,
                                                     int taxon,
                                                     GradientComponents components,
                                                     double[] transposedLoadings,
                                                     double[] rawGamma,
                                                     double[][] gradArray,
                                                     int offset) {
        double[] fty = components.fty;
        double[] ftfl = components.ftfl;


        for (int trait = 0; trait < dimTrait; ++trait) {
            int dataInd = taxon * dimTrait + trait;
            if (factorAnalysisLikelihood.getDataMissingIndicators()[dataInd]) {
                continue;
            }
            double dat = data[dataInd];
            gradArray[index][offset + trait] += 0.5 * (1 / rawGamma[trait] - dat * dat);

            for (int factor = 0; factor < dimFactors; ++factor) {
                int loadingsInd = trait * dimFactors + factor;
                int ind = factor * dimTrait + trait;

                gradArray[index][offset + trait] +=
                        (fty[ind] - 0.5 * ftfl[ind]) * transposedLoadings[loadingsInd];
            }
//            gradArray[index][offset + trait] +=
//                    (fty[ind] - 0.5 * ftfl[ind]) * transposedLoadings[ind];
        }
//        for (int factor = 0; factor < dimFactors; ++factor) {
//            for (int trait = 0; trait < dimTrait; ++trait) {
//                int ind = factor * dimTrait + trait;
//                gradArray[index][offset + trait] +=
//                        (fty[ind] - 0.5 * ftfl[ind]) * transposedLoadings[ind];
//            }
//        }
//        for (int trait = 0; trait < dimTrait; ++trait) {
//            double dat = data[taxon * dimTrait + trait];
//            gradArray[index][offset + trait] += 0.5 * (1 / rawGamma[trait] - dat * dat); //TODO: need to deal w/ missing data
//        }
    }

    @Override
    protected void computeGradientForOneTaxon(final int index,
                                              final int taxon,
                                              final ReadableMatrix loadings,
                                              final double[] transposedLoadings,
                                              final ReadableVector gamma,
                                              final double[] rawGamma,
                                              final WrappedNormalSufficientStatistics statistic,
                                              final double[][] gradArray) {

        if (TIMING) {
            stopWatches[0].start();
        }

        final MeanAndMoment meanAndMoment = getMeanAndMoment(taxon, statistic);

        if (TIMING) {
            stopWatches[0].stop();
            stopWatches[1].start();
        }

        GradientComponents components = computeGradientComponents(taxon, transposedLoadings, meanAndMoment);
        computeLoadingsGradientForOneTaxon(index, components, rawGamma, gradArray);
        computePrecisionGradientForOneTaxon(index, taxon, components, transposedLoadings, rawGamma, gradArray, dimFactors * dimTrait);


        if (TIMING) {
            stopWatches[1].stop();
        }
//        }
    }


}
