package dr.evomodel.coalescent;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.Reportable;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * @author Xiang Ji
 */

public class MultilocusNPCoalescentLikelihoodGradient extends AbstractModel implements
        GradientWrtParameterProvider, Reportable {

    private final MultilocusNonparametricCoalescentLikelihood likelihood;
    private final Parameter parameter;
    private final GradientProvider provider;
    private boolean gradientKnown;
    private double[] gradient;

    public MultilocusNPCoalescentLikelihoodGradient(MultilocusNonparametricCoalescentLikelihood likelihood,
                                                    Parameter parameter) {
        super("MultilocusNPCoalescentLikelihoodGradient");

        this.likelihood = likelihood;
        this.parameter = parameter;

        addVariable(parameter);
        addModel(likelihood);

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
        this.gradientKnown = false;
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
        if(!gradientKnown) {
            gradient = provider.getGradientLogDensity(likelihood.getLogPopSizes().getParameterValues());
            gradientKnown = true;
        }
        return gradient;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof MultilocusNonparametricCoalescentLikelihood) {
            gradientKnown = false;
        } else {
            throw new RuntimeException("Unknown object");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        gradientKnown = false;
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }
}
