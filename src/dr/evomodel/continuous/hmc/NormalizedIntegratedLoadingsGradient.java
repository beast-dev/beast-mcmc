package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledMatrixParameter;
import dr.util.TaskPool;

public class NormalizedIntegratedLoadingsGradient extends IntegratedLoadingsGradient {

    private final ScaledMatrixParameter scaledMatrix;

    public NormalizedIntegratedLoadingsGradient(TreeDataLikelihood treeDataLikelihood,
                                                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                                TaskPool taskPool,
                                                ThreadUseProvider threadUseProvider,
                                                RemainderCompProvider remainderCompProvider) {
        super(treeDataLikelihood, likelihoodDelegate, factorAnalysisLikelihood, taskPool, threadUseProvider,
                remainderCompProvider);

        this.scaledMatrix = (ScaledMatrixParameter) factorAnalysisLikelihood.getLoadings();

    }

    @Override
    public Parameter getParameter() {
        return scaledMatrix.getMatrixParameter();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] gradLogDensity = super.getGradientLogDensity();
        int offset = 0;
        for (int factor = 0; factor < dimFactors; factor++) {
            double scale = scaledMatrix.getScaleParameter().getParameterValue(factor);

            for (int trait = 0; trait < dimTrait; trait++) {
                gradLogDensity[offset + trait] *= scale;
            }
            offset += dimTrait;
        }

        return gradLogDensity;
    }

}
