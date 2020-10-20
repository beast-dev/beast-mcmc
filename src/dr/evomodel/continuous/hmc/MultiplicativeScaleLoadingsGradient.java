package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledMatrixParameter;
import dr.util.TaskPool;


public class MultiplicativeScaleLoadingsGradient extends ScaleIntegratedLoadingsGradient {


    private final Parameter scale;
    private final Parameter multipliers;

    public MultiplicativeScaleLoadingsGradient(TreeDataLikelihood treeDataLikelihood,
                                               ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                               IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                               TaskPool taskPool,
                                               ThreadUseProvider threadUseProvider,
                                               RemainderCompProvider remainderCompProvider,
                                               Parameter parameter) {
        super(treeDataLikelihood, likelihoodDelegate, factorAnalysisLikelihood, taskPool, threadUseProvider, remainderCompProvider);
        this.scale = ((ScaledMatrixParameter) factorAnalysisLikelihood.getLoadings()).getScaleParameter();
        this.multipliers = parameter;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] scaledGradient = super.getGradientLogDensity();
        double[] transformedGradient = new double[scaledGradient.length];

        for (int i = 0; i < scaledGradient.length; i++) {
            for (int j = i; j < scaledGradient.length; j++) {
                transformedGradient[i] += scaledGradient[j] * scale.getParameterValue(j) / multipliers.getParameterValue(i);
            }
        }

        return transformedGradient;
    }

    @Override
    public Parameter getParameter() {
        return multipliers;
    }

}
