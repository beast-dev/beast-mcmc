/*
 * MultivariateTDiffusionModel.java
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

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.distributions.TDistribution;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class MultivariateTDiffusionModel extends MultivariateDiffusionModel {

    public static final String DIFFUSION_PROCESS = "multivariateTDiffusionModel";
    public static final String DIFFUSION_CONSTANT = "precisionParameter";
    //	public static final String BIAS = "mu";
    //	public static final String PRECISION_TREE_ATTRIBUTE = "precision";
    public static final String DF = "dfParameter";

    public MultivariateTDiffusionModel(Parameter df, MatrixParameterInterface precision) {
        super();
        this.dfParameter = df;
        this.precisionParameter = precision;
        addVariable(dfParameter);
        addVariable(precisionParameter);

    }

    public double getLogLikelihood(double[] start, double[] stop, double time) {

        double df = dfParameter.getParameterValue(0);
//		double scale0 = diffusionPrecisionMatrixParameter.getParameterValue(0,0);
//		double scale1 = diffusionPrecisionMatrixParameter.getParameterValue(1,1);
        double scale0 = precisionParameter.getParameterValue(0);
        double scale1 = precisionParameter.getParameterValue(1);
        // todo Make this truely multivariate
        return TDistribution.logPDF(stop[0], start[0], time / scale0, df) +
                TDistribution.logPDF(stop[1], start[1], time / scale1, df);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    // No need to just call super of function, is there?

//    public void handleModelChangedEvent(Model model, Object object, int index) {
//        super.handleModelChangedEvent(model, object, index);    //To change body of overridden methods use File | Settings | File Templates.
//    }
/*
    protected final void handleVariableChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
//		if( parameter.getId().compareTo("allTraits")!= 0)
//			System.err.println("parameter = "+parameter.getId());
        super.handleVariableChangedEvent(parameter, index, type);    //To change body of overridden methods use File | Settings | File Templates.
    }
*/

    public MatrixParameterInterface getPrecisionParameter() {
        return precisionParameter;
    }

    protected void calculatePrecisionInfo() {
        // Nothing to do?
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(DIFFUSION_CONSTANT);
            MatrixParameterInterface diffusionParam = (MatrixParameterInterface) cxo.getChild(MatrixParameterInterface.class);
            cxo = xo.getChild(DF);
            Parameter df = (Parameter) cxo.getChild(Parameter.class);

            return new MultivariateTDiffusionModel(df, diffusionParam);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Describes a multivariate t-distribution diffusion process.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(DIFFUSION_CONSTANT,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(DF,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
        };

        public Class getReturnType() {
            return MultivariateTDiffusionModel.class;
        }
    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Parameter dfParameter;
    private final MatrixParameterInterface precisionParameter;

}
