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
    private final SubstitutionModel substitutionModel;

    private final WrtParameter wrtParameter;
    // private final AbstractGlmSubstitutionModelGradient.ParameterMap parameterMap;


    public StructuredCoalescentLikelihoodGradient(BastaLikelihood BastaLikelihood,
                                                  SubstitutionModel substitutionModel,
                                                  WrtParameter wrtParameter) {
//        this.structuredCoalescentLikelihood = structuredCoalescentLikelihood;
        this.structuredCoalescentLikelihood = BastaLikelihood;
        this.substitutionModel = substitutionModel;
        // this.parameterMap = makeParameterMap(glm);
        this.wrtParameter = wrtParameter;

    }



    @Override
    public Likelihood getLikelihood() {
        return structuredCoalescentLikelihood;
    }

    @Override
    public Parameter getParameter() {
            return wrtParameter.getParameter(structuredCoalescentLikelihood, substitutionModel);
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return wrtParameter.getGradientLogDensity(structuredCoalescentLikelihood);
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

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1E-1);
        sb.append(message);


        return  sb.toString();
    }

    public enum WrtParameter {
        MIGRATION_RATE("migrationRate") {
            @Override
            double[] getGradientLogDensity(BastaLikelihood structuredCoalescentLikelihood) {
                return structuredCoalescentLikelihood.getGradientLogDensity();
//            double[] getGradientLogDensity(StructuredCoalescentLikelihood BastaLikelihood) {
//                return BastaLikelihood.getGradientLogDensity();
            }

            @Override
            Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                assert(substitutionModel instanceof SVSComplexSubstitutionModel);
                SVSComplexSubstitutionModel svsComplexSubstitutionModel = (SVSComplexSubstitutionModel) substitutionModel;
                return svsComplexSubstitutionModel.getRatesParameter();
            }

        },

        POPULATION_SIZE("populationSize") {
            @Override
            double[] getGradientLogDensity(BastaLikelihood structuredCoalescentLikelihood) {
                return structuredCoalescentLikelihood.getPopSizeGradientLogDensity();
            }
            

            @Override
            Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                return structuredCoalescentLikelihood.getPopSizes();
            }
        };

        WrtParameter(String name) {
            this.name = name;
        }

        abstract double[] getGradientLogDensity(BastaLikelihood structuredCoalescentLikelihood);

        abstract Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel);

        private final String name;

        public static WrtParameter factory(String match) {
            for (WrtParameter type : WrtParameter.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }
    }

}
