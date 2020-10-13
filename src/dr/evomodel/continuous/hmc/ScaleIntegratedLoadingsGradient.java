package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledMatrixParameter;
import dr.util.TaskPool;

public class ScaleIntegratedLoadingsGradient extends IntegratedLoadingsGradient {

    private final ScaledMatrixParameter scaledMatrix;

    public ScaleIntegratedLoadingsGradient(TreeDataLikelihood treeDataLikelihood,
                                           ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                           IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                           TaskPool taskPool,
                                           IntegratedLoadingsGradient.ThreadUseProvider threadUseProvider,
                                           IntegratedLoadingsGradient.RemainderCompProvider remainderCompProvider) {
        super(treeDataLikelihood, likelihoodDelegate, factorAnalysisLikelihood, taskPool, threadUseProvider,
                remainderCompProvider);

        this.scaledMatrix = (ScaledMatrixParameter) factorAnalysisLikelihood.getLoadings();

    }

    @Override
    public Parameter getParameter() {
        return scaledMatrix.getScaleParameter();
    }

    @Override
    public int getDimension() {
        return dimFactors;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] loadingsGradient = null;
        try {
            loadingsGradient = super.getGradientLogDensity();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }

        double[] scaleGradient = new double[dimFactors];

        int offset = 0;
        for (int factor = 0; factor < dimFactors; factor++) {

            for (int trait = 0; trait < dimTrait; trait++) {
                scaleGradient[factor] += loadingsGradient[offset + trait] *
                        scaledMatrix.getMatrixParameter().getParameterValue(trait, factor);
            }
            offset += dimTrait;
        }

        return scaleGradient;
    }
}
