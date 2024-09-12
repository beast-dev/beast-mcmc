/*
 * ParameterValuesReport.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

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
