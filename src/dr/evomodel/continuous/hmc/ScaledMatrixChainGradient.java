package dr.evomodel.continuous.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledMatrixParameter;
import dr.xml.*;

public class ScaledMatrixChainGradient implements GradientWrtParameterProvider {

    private final GradientWrtParameterProvider originalGradient;
    private final ScaledMatrixParameter parameter;
    private final ComponentProvider componentProvider;

    ScaledMatrixChainGradient(GradientWrtParameterProvider originalGradient, ComponentProvider componentProvider) {
        this.originalGradient = originalGradient;
        this.parameter = (ScaledMatrixParameter) originalGradient.getParameter();
        this.componentProvider = componentProvider;
    }

    @Override
    public Likelihood getLikelihood() {
        return originalGradient.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return componentProvider.getParameter(parameter);
    }

    @Override
    public int getDimension() {
        return componentProvider.getDimension(parameter);
    }

    @Override
    public double[] getGradientLogDensity() {
        return componentProvider.chainGradient(originalGradient.getGradientLogDensity(), parameter);
    }

    public enum ComponentProvider {
        MATRIX("matrix") {
            @Override
            Parameter getParameter(ScaledMatrixParameter parameter) {
                return parameter.getMatrixParameter();
            }

            @Override
            double[] chainGradient(double[] gradient, ScaledMatrixParameter parameter) {
                int offset = 0;
                int nCols = parameter.getRowDimension();
                int nRows = parameter.getColumnDimension();
                for (int factor = 0; factor < nRows; factor++) {
                    double scale = parameter.getScaleParameter().getParameterValue(factor);

                    for (int trait = 0; trait < nCols; trait++) {
                        gradient[offset + trait] *= scale;
                    }
                    offset += nCols;
                }

                return gradient;
            }
        },

        SCALE("scale") {
            @Override
            Parameter getParameter(ScaledMatrixParameter parameter) {
                return parameter.getScaleParameter();
            }

            @Override
            double[] chainGradient(double[] gradient, ScaledMatrixParameter parameter) {
                int nCols = parameter.getRowDimension();
                int nRows = parameter.getColumnDimension();

                double[] scaleGradient = new double[nRows];

                int offset = 0;
                for (int factor = 0; factor < nRows; factor++) {

                    for (int trait = 0; trait < nCols; trait++) {
                        scaleGradient[factor] += gradient[offset + trait] *
                                parameter.getMatrixParameter().getParameterValue(trait, factor);
                    }
                    offset += nCols;
                }

                return scaleGradient;
            }
        };

        public final String name;

        ComponentProvider(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(ScaledMatrixParameter parameter);

        abstract double[] chainGradient(double[] gradient, ScaledMatrixParameter parameter);


        public int getDimension(ScaledMatrixParameter parameter) {
            return getParameter(parameter).getDimension();
        }


    }

    private static final String SCALED_GRADIENT = "scaledMatrixGradient";
    private static final String COMPONENT = "component";
    private static final String SCALE = "scale";
    private static final String MATRIX = "matrix";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String component = xo.getStringAttribute(COMPONENT);
            ComponentProvider provider = null;
            for (ComponentProvider componentProvider : ComponentProvider.values()) {
                if (component.equalsIgnoreCase(componentProvider.name)) {
                    provider = componentProvider;
                }
            }
            if (provider == null) {
                throw new XMLParseException("Unrecognized 'component'. Must be '" + SCALE + "' or '" + MATRIX + "'.");
            }

            GradientWrtParameterProvider gradient =
                    (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
            return new ScaledMatrixChainGradient(gradient, provider);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(GradientWrtParameterProvider.class),
                    AttributeRule.newStringRule(COMPONENT)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return ScaledMatrixChainGradient.class;
        }

        @Override
        public String getParserName() {
            return SCALED_GRADIENT;
        }
    };
}
