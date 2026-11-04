package dr.evomodel.continuous;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.xml.Reportable;

public class SampledCategoricalVarianceGradient extends AbstractModel
        implements GradientWrtParameterProvider, Reportable {

    private final SampledCategoricalVarianceModel likelihood;
    private final Parameter parameter;
    private final GradientWrtParameterProvider gradientProvider;

    public SampledCategoricalVarianceGradient(String name,
                                              SampledCategoricalVarianceModel likelihood,
                                              Parameter parameter) {
        super(name);

        this.likelihood = likelihood;
        this.parameter = parameter;

        addModel(likelihood);
        addVariable(parameter);

        if (parameter != likelihood.getTraitsParameter()) {
            throw new IllegalArgumentException("Unknown parameter");
        }

        this.gradientProvider = makeGradientProviderWrtTraits();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    GradientWrtParameterProvider makeGradientProviderWrtTraits() {

        final Parameter traits = likelihood.getTraitsParameter();
        final Parameter means = likelihood.getMeansParameter();
        final SampledCategoricalVarianceModel.Parametrization parametrization = likelihood.getParameterization();

        return new GradientWrtParameterProvider() {
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
                return getParameter().getDimension();
            }

            @Override
            public double[] getGradientLogDensity() {

                double[] gradient = new double[getDimension()];
                for (int i = 0; i < getDimension(); ++i) {
                    double sd = parametrization.getSd(i);
                    double mean = means.getParameterValue(i);
                    if (Double.isFinite(mean)) {
                        gradient[i] = NormalDistribution.gradLogPdf(traits.getParameterValue(i),
                                mean, sd);
                    } else {
                        gradient[i] = 0.0;
                    }
                }

                return gradient;
            }
        };
    }

    @Override
    public Likelihood getLikelihood() {
        return gradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return gradientProvider.getParameter();
    }

    @Override
    public int getDimension() {
        return gradientProvider.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return gradientProvider.getGradientLogDensity();
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                getParameter().getBounds().getLowerLimit(0),
                getParameter().getBounds().getUpperLimit(0), null);
    }
}
