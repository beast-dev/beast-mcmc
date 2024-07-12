package dr.evomodel.coalescent.basta;

import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.xml.Reportable;

public class StructuredCoalescentLikelihoodGradient implements
        GradientWrtParameterProvider, ModelListener, Reportable, Loggable {

    private final BastaLikelihood structuredCoalescentLikelihood;
    private final WrtParameter wrtParameter;

    private final Parameter parameter;
    private final Parameter chainRuleDependent;

    private final int stateCount;

    public StructuredCoalescentLikelihoodGradient(BastaLikelihood BastaLikelihood,
                                                  SubstitutionModel substitutionModel,
                                                  WrtParameter wrtParameter) {
        this.structuredCoalescentLikelihood = BastaLikelihood;
        this.wrtParameter = wrtParameter;

        this.parameter = wrtParameter.getParameter(structuredCoalescentLikelihood, substitutionModel);
        this.chainRuleDependent = wrtParameter.getChainRuleDependent(structuredCoalescentLikelihood, substitutionModel);

        this.stateCount = structuredCoalescentLikelihood.getSubstitutionModel().getFrequencyModel().getFrequencyCount();
    }

    @Override
    public Likelihood getLikelihood() {
        return structuredCoalescentLikelihood;
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
        return structuredCoalescentLikelihood.getGradientLogDensity(this);
    }

    double[] chainRule(double[] gradient) {
        return wrtParameter.chainRule(gradient, chainRuleDependent);
    }

    boolean requiresTransitionMatrices() {
        return wrtParameter.requiresTransitionMatrices();
    }

    WrtParameter getType() { return wrtParameter; }

    public int getIntermediateGradientDimension() {
//        return structuredCoalescentLikelihood.
        return wrtParameter.getIntermediateGradientDimension(stateCount);
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    public void modelRestored(Model model) {

    }


    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 10.0);
        sb.append(message);


        return  sb.toString();
    }

    public enum WrtParameter {
        MIGRATION_RATE("migrationRate") {
            @Override
            Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                assert(substitutionModel instanceof SVSComplexSubstitutionModel);
                SVSComplexSubstitutionModel svsComplexSubstitutionModel = (SVSComplexSubstitutionModel) substitutionModel;
                return svsComplexSubstitutionModel.getRatesParameter();
            }

            @Override
            double[] chainRule(double[] gradient, Parameter parameter) {
                final int dim = parameter.getDimension();

                double[] chainedGradient = new double[dim * (dim - 1)];

                int k = 0;
                for (int i = 0; i < dim; ++i) {
                    for (int j = i + 1; j < dim; ++j) {
                        chainedGradient[k] = (gradient[i * dim + j] - gradient[i * dim + i]) * parameter.getParameterValue(j) ;
                        ++k;
                    }
                }

                for (int j = 0; j < dim; ++j) {
                    for (int i = j + 1; i < dim; ++i) {
                        chainedGradient[k] =(gradient[i * dim + j] - gradient[i * dim + i]) * parameter.getParameterValue(j);
                        ++k;
                    }
                }
                return chainedGradient;
            }

            @Override
            int getIntermediateGradientDimension(int stateCount) {
                return stateCount * stateCount;
            }

            @Override
            boolean requiresTransitionMatrices() {
                return false;
            }

            @Override
            public Parameter getChainRuleDependent(BastaLikelihood structuredCoalescentLikelihood,
                                                   SubstitutionModel substitutionModel) {
                return substitutionModel.getFrequencyModel().getFrequencyParameter();
            }

        },

        POPULATION_SIZE("populationSize") {
            @Override
            Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                return structuredCoalescentLikelihood.getPopSizes();
            }

            @Override
            double[] chainRule(double[] gradient, Parameter parameter) {
                final int dim = parameter.getDimension();

                for (int i = 0; i < dim; ++i) {//
                    double popSize = parameter.getParameterValue(i);
                    gradient[i] /= -(popSize * popSize);
                }

                return gradient;
            }

            @Override
            int getIntermediateGradientDimension(int stateCount) {
                return stateCount;
            }

            @Override
        boolean requiresTransitionMatrices() {
            return false;
        }

            @Override
            public Parameter getChainRuleDependent(BastaLikelihood structuredCoalescentLikelihood,
                                                   SubstitutionModel substitutionModel) {
                return structuredCoalescentLikelihood.getPopSizes();
            }
        };

        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel);

        abstract double[] chainRule(double[] gradient, Parameter parameter);

        abstract int getIntermediateGradientDimension(int stateCount);

        abstract boolean requiresTransitionMatrices();

        private final String name;

        public static WrtParameter factory(String match) {
            for (WrtParameter type : WrtParameter.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }

        public abstract Parameter getChainRuleDependent(BastaLikelihood structuredCoalescentLikelihood,
                                                        SubstitutionModel substitutionModel);
    }
}
