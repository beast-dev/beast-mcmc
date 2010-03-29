/*
 * SetOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
 */

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A generic operator for selecting uniformly from a discrete set of values.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: SetOperator.java,v 1.12 2005/05/24 20:26:00 rambaut Exp $
 */
public class SetOperator extends SimpleMCMCOperator {

    public SetOperator(Parameter parameter, double[] values) {
        this.parameter = parameter;
        this.values = values;
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        int index = MathUtils.nextInt(values.length);
        double newValue = values[index];

        if (newValue < parameter.getBounds().getLowerLimit(index) || newValue > parameter.getBounds().getUpperLimit(index)) {
            throw new OperatorFailedException("proposed value outside boundaries");
        }

        parameter.setParameterValue(index, newValue);

        return 0.0;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return "setOperator";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double[] values = xo.getDoubleArrayAttribute("set");
            double weight = xo.getDoubleAttribute("weight");

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            System.out.println("Creating set operator for parameter " + parameter.getParameterName());
            System.out.print("  set = {" + values[0]);
            for (int i = 1; i < values.length; i++) {
                System.out.print(", " + values[i]);
            }
            System.out.println("}");

            SetOperator operator = new SetOperator(parameter, values);
            operator.setWeight(weight);

            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an operator on a set.";
        }

        public Class getReturnType() {
            return SetOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleArrayRule("set"),
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(Parameter.class)
        };
    };

    public Element createOperatorElement(Document document) {
        throw new RuntimeException("Not implememented!");
    }

    public String getOperatorName() {
        return "setOperator(" + parameter.getParameterName() + ")";
    }

    public String getPerformanceSuggestion() {
        return "No suggestions";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private double[] values;

}
