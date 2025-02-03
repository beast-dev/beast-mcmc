package dr.evomodel.coalescent;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * @author Xiang Ji
 */

public class MultilocusNPCoalescentLikelihoodGradient implements GradientWrtParameterProvider, Reportable {

    private final MultilocusNonparametricCoalescentLikelihood likelihood;
    private final Parameter parameter;
    private final GradientProvider provider;

    public MultilocusNPCoalescentLikelihoodGradient(MultilocusNonparametricCoalescentLikelihood likelihood,
                               Parameter parameter) {

        this.likelihood = likelihood;
        this.parameter = parameter;

        if (parameter == likelihood.getLogPopSizes()) {
            provider = new GradientProvider() {
                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return likelihood.getGradientLogDensity(x);
                }
            };
        }
        else {
            throw new IllegalArgumentException("Not yet implemented");
//            provider = likelihood.getGradientWrt(parameter);
        }
    }
    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, null);
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        // TODO Can cache here
        return provider.getGradientLogDensity(likelihood.getLogPopSizes().getParameterValues());
    }
}
