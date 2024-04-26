package dr.evomodel.treedatalikelihood.discrete;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.MaskedGradient;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Likelihood;
import dr.inference.model.MaskedParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.Set;

public class MaskedRatioGradient implements GradientWrtParameterProvider {

    final private GradientWrtParameterProvider nodeHeightGradient;
    final private NodeHeightTransform transform;
    final private Parameter nodeHeightParameter;
    final private WeightedLikelihood weightedLikelihood;
    final private Parameter ratios;
    final private Parameter mask;
    final private MaskedGradient maskedGradient;
    final private MaskedParameter maskedParameter;

    public MaskedRatioGradient(NodeHeightTransform transform,
                               GradientWrtParameterProvider nodeHeightGradient,
                               Parameter ratios,
                               Parameter nodeHeightParameter,
                               Parameter mask) {

        this.nodeHeightGradient = nodeHeightGradient;
        this.transform = transform;
        this.nodeHeightParameter = nodeHeightParameter;
        this.mask = mask;
        this.ratios = ratios;
        this.weightedLikelihood = new WeightedLikelihood();
        this.maskedGradient = new MaskedGradient(nodeHeightGradient, mask);
        this.maskedParameter = new MaskedParameter(ratios, mask, true);
    }

    @Override
    public Likelihood getLikelihood() {
        return weightedLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return nodeHeightParameter;
    }

    @Override
    public int getDimension() {
        return maskedParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return maskedGradient.getGradientLogDensity();
    }

    private class WeightedLikelihood implements Likelihood {

        @Override
        public LogColumn[] getColumns() {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public Model getModel() {
            throw new RuntimeException("not yet implemented");
        }

        @Override
        public double getLogLikelihood() {
            return nodeHeightGradient.getLikelihood().getLogLikelihood() + transform.getLogJacobian(nodeHeightParameter.getParameterValues());
        }

        @Override
        public void makeDirty() {

        }

        @Override
        public String prettyName() {
            return null;
        }

        @Override
        public Set<Likelihood> getLikelihoodSet() {
            return null;
        }

        @Override
        public boolean isUsed() {
            return false;
        }

        @Override
        public void setUsed() {

        }

        @Override
        public boolean evaluateEarly() {
            return false;
        }

        @Override
        public String getId() {
            return null;
        }

        @Override
        public void setId(String id) {

        }
    }
}
