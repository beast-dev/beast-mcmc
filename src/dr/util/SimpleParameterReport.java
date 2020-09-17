package dr.util;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.*;

public class SimpleParameterReport implements Reportable {

    private final Parameter parameter;

    SimpleParameterReport(Parameter parameter) {
        this.parameter = parameter;
    }


    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder(parameter.getVariableName() + " values:");
        if (parameter instanceof MatrixParameterInterface) {
            for (int i = 0; i < ((MatrixParameterInterface) parameter).getRowDimension(); i++) {
                sb.append("\n\t");
                for (int j = 0; j < ((MatrixParameterInterface) parameter).getColumnDimension(); j++) {
                    sb.append(((MatrixParameterInterface) parameter).getParameterValue(i, j));
                    sb.append(" ");
                }
            }
        } else {
            for (int i = 0; i < parameter.getDimension(); i++) {
                sb.append(" ");
                sb.append(parameter.getParameterValue(i));
            }
        }
        sb.append("\n\n");
        return sb.toString();
    }


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String SIMPLE_PARAMETER_REPORT = "simpleParameterReport";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new SimpleParameterReport(parameter);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Report with only the values of a parameter.";
        }

        @Override
        public Class getReturnType() {
            return SimpleParameterReport.class;
        }

        @Override
        public String getParserName() {
            return SIMPLE_PARAMETER_REPORT;
        }
    };

}


