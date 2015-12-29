/*
 * BiasedMultivariateDiffusionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.xml.*;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evolution.tree.Tree;

/**
 * @author Marc A. Suchard
 */

public class BiasedMultivariateDiffusionModel extends MultivariateDiffusionModel {

    public static final String BIASED_DIFFUSION_PROCESS = "biasedMultivariateDiffusionModel";
    public static final String BIAS_PARAMETER = "biasParameter";
    public static final String BIAS_TREE_ATTRIBUTE = "bias";

    BiasedMultivariateDiffusionModel(Parameter biasParam, MatrixParameter diffusionParam) {
        super(diffusionParam);
        this.biasParam = biasParam;
        addVariable(biasParam);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == diffusionPrecisionMatrixParameter)
            calculatePrecisionInfo();
        // else is bias and do nothing
    }

    @Override
    protected double calculateLogDensity(double[] start, double[] stop, double time) {

        final int dim = stop.length;
        double[] bias = biasParam.getParameterValues();
        for (int i = 0; i < dim; i++) {
            bias[i] *= time;
            bias[i] += start[i];
        }



        return super.calculateLogDensity(bias, stop, time);
    }

    public String[] getTreeAttributeLabel() {
        return new String[]{
                PRECISION_TREE_ATTRIBUTE,
                BIAS_TREE_ATTRIBUTE
        };
    }

    public String[] getAttributeForTree(Tree tree) {
        return new String[]{
                diffusionPrecisionMatrixParameter.toSymmetricString(),
                toParameterString(biasParam)
        };
    }

    private String toParameterString(Parameter param) {
        StringBuffer sb = new StringBuffer("{");
        int dim = param.getDimension();
        for (int i = 0; i < dim; i++) {
            sb.append(String.format("%5.4e", param.getParameterValue(i)));
            if (i < dim - 1)
                sb.append(",");
        }
        sb.append("}");
        return sb.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BIASED_DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(BIAS_PARAMETER);
            Parameter biasParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(DIFFUSION_CONSTANT);
            MatrixParameter diffusionParam = (MatrixParameter) cxo.getChild(MatrixParameter.class);

            return new BiasedMultivariateDiffusionModel(biasParam, diffusionParam);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Describes a multivariate normal diffusion process.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(DIFFUSION_CONSTANT,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
                new ElementRule(BIAS_PARAMETER,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        };

        public Class getReturnType() {
            return MultivariateDiffusionModel.class;
        }
    };

    private Parameter biasParam;

}
