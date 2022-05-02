package dr.xml.unittest;

import dr.inference.model.Parameter;
import dr.xml.*;

public class ParameterValuesReport implements Reportable {

    private final Parameter parameter;

    public ParameterValuesReport(Parameter parameter) {
        this.parameter = parameter;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < parameter.getDimension(); i++) {
            sb.append(parameter.getParameterValue(i));
            sb.append(" ");
        }
        return sb.toString();
    }

    private static final String PARAMETER_VALUES = "parameterValues";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter param = (Parameter) xo.getChild(Parameter.class);
            return new ParameterValuesReport(param);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "prints parameter values to string";
        }

        @Override
        public Class getReturnType() {
            return ParameterValuesReport.class;
        }

        @Override
        public String getParserName() {
            return PARAMETER_VALUES;
        }
    };
}
