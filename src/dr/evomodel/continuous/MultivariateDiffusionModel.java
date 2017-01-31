/*
 * MultivariateDiffusionModel.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.inference.model.*;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc Suchard
 */


public class MultivariateDiffusionModel extends AbstractModel implements TreeAttributeProvider {

    public static final String DIFFUSION_PROCESS = "multivariateDiffusionModel";
    public static final String DIFFUSION_CONSTANT = "precisionMatrix";
    public static final String PRECISION_TREE_ATTRIBUTE = "precision";

    public static final double LOG2PI = Math.log(2*Math.PI);

    /**
     * Construct a diffusion model.
     */

    public MultivariateDiffusionModel(MatrixParameterInterface diffusionPrecisionMatrixParameter) {

        super(DIFFUSION_PROCESS);

        this.diffusionPrecisionMatrixParameter = diffusionPrecisionMatrixParameter;
        calculatePrecisionInfo();
        addVariable(diffusionPrecisionMatrixParameter);

    }

    public MultivariateDiffusionModel() {
        super(DIFFUSION_PROCESS);
    }


//    public void randomize(Parameter trait) {
//    }

    public void check(Parameter trait) throws XMLParseException {
        assert trait != null;
    }

    public MatrixParameterInterface getPrecisionParameter() {
        checkVariableChanged();
        return diffusionPrecisionMatrixParameter;
    }

    public double[][] getPrecisionmatrix() {
        if (diffusionPrecisionMatrixParameter != null) {
            checkVariableChanged();
            return diffusionPrecisionMatrixParameter.getParameterAsMatrix();
        }
        return null;
    }

    public double getDeterminantPrecisionMatrix() {
        checkVariableChanged();
        return determinatePrecisionMatrix;
    }

    /**
     * @return the log likelihood of going from start to stop in the given time
     */
    public double getLogLikelihood(double[] start, double[] stop, double time) {

        if (time == 0) {
            boolean equal = true;
            for(int i=0; i<start.length; i++) {
                if( start[i] != stop[i] ) {
                    equal = false;
                    break;
                }
            }
            if (equal)
                return 0.0;
            return Double.NEGATIVE_INFINITY;
        }

        return calculateLogDensity(start, stop, time);
    }

    protected void checkVariableChanged(){
        if(variableChanged){
            calculatePrecisionInfo();
            variableChanged=false;
        }
    }

    protected double calculateLogDensity(double[] start, double[] stop, double time) {
        checkVariableChanged();
        final double logDet = Math.log(determinatePrecisionMatrix);
        return MultivariateNormalDistribution.logPdf(stop, start, diffusionPrecisionMatrix, logDet, time);
    }

    // todo should be a test, no?
    public static void main(String[] args) {
        double[] start = {1, 2};
        double[] stop = {0, 0};
        double[][] precision = {{2, 0.5}, {0.5, 1}};
        double scale = 0.2;
        MatrixParameter precMatrix = new MatrixParameter("Hello");
        precMatrix.addParameter(new Parameter.Default(precision[0]));
        precMatrix.addParameter(new Parameter.Default(precision[1]));
        MultivariateDiffusionModel model = new MultivariateDiffusionModel(precMatrix);
        System.err.println("logPDF = " + model.calculateLogDensity(start, stop, scale));
        System.err.println("Should be -19.948");
    }

    protected void calculatePrecisionInfo() {
        diffusionPrecisionMatrix = diffusionPrecisionMatrixParameter.getParameterAsMatrix();
        determinatePrecisionMatrix =
                MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(
                        diffusionPrecisionMatrix);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************
    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        variableChanged=true;
//        calculatePrecisionInfo();
    }

    protected void storeState() {
        savedDeterminatePrecisionMatrix = determinatePrecisionMatrix;
        savedDiffusionPrecisionMatrix = diffusionPrecisionMatrix;
        storedVariableChanged=variableChanged;
    }

    protected void restoreState() {
        determinatePrecisionMatrix = savedDeterminatePrecisionMatrix;
        diffusionPrecisionMatrix = savedDiffusionPrecisionMatrix;
        variableChanged=storedVariableChanged;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    public String[] getTreeAttributeLabel() {
        return new String[] {PRECISION_TREE_ATTRIBUTE};
    }

    public String[] getAttributeForTree(Tree tree) {
        if (diffusionPrecisionMatrixParameter != null) {
            return new String[] {diffusionPrecisionMatrixParameter.toSymmetricString()};
        }

        diffusionPrecisionMatrixParameter.toString();
        return new String[] { "null" };
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(DIFFUSION_CONSTANT);
            MatrixParameterInterface diffusionParam = (MatrixParameterInterface)
                    cxo.getChild(MatrixParameterInterface.class);

            return new MultivariateDiffusionModel(diffusionParam);
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
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
        };

        public Class getReturnType() {
            return MultivariateDiffusionModel.class;
        }
    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    protected MatrixParameterInterface diffusionPrecisionMatrixParameter;
    private double determinatePrecisionMatrix;
    private double savedDeterminatePrecisionMatrix;
    private double[][] diffusionPrecisionMatrix;
    private double[][] savedDiffusionPrecisionMatrix;

    private boolean variableChanged=true;
    private boolean storedVariableChanged;

}

