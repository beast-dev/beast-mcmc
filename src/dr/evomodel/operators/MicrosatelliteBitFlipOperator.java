/*
 * MsatBitFlipOperator.java
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

package dr.evomodel.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 *  This operator performs bitflip operation on the bit vector representing the model.
 * 
 */
public class MicrosatelliteBitFlipOperator extends SimpleMCMCOperator {
    private Parameter parameter;
    private Parameter dependencies;
    private int[] variableIndices;
    public static final int PRESENT = 1;
    public static final int ABSENT = 0;
    public static final int NO_DEPENDENCY = -1;
    public static final String MODEL_CHOOSE = "modelChoose";
    public static final String DEPENDENCIES = "dependencies";
    public static final String VARIABLE_INDICES = "variableIndices";


    public MicrosatelliteBitFlipOperator(Parameter parameter, Parameter dependencies, double weight, int[] variableIndices){
        this.parameter = parameter;
        this.dependencies = dependencies;
        this.variableIndices = variableIndices;
        if(parameter.getDimension() != dependencies.getDimension())
            throw new RuntimeException("Dimenension of the parameter ("+parameter.getDimension()+
                    ") does not equal to the dimension of the dependencies parameter("+dependencies.getDimension()+").");
        setWeight(weight);
    }

    public String getOperatorName(){
        return "msatModelSwitch(" + parameter.getParameterName() + ")";
    }

    public double doOperation() {

        double logq = 0.0;
        double[] bitVec = new double[parameter.getDimension()];
        for(int i = 0; i < bitVec.length; i++){
            bitVec[i] = parameter.getParameterValue(i);
        }
        //int index = (int)MathUtils.nextDouble()*parameter.getDimension();
        int index = variableIndices[MathUtils.nextInt(variableIndices.length)];
        //System.out.println(index);
        int oldVal  = (int)parameter.getParameterValue(index);
        int newVal = -1;
        if(oldVal == ABSENT){
            newVal = PRESENT;
        }else if(oldVal == PRESENT){
           newVal = ABSENT;
        }else{
            throw new RuntimeException("The parameter can only take values 0 or 1.");
        }
        bitVec[index] = newVal;
        for(int i = 0; i < bitVec.length; i++){
            int dependentInd = (int)dependencies.getParameterValue(i);
            if(dependentInd > NO_DEPENDENCY){
                if(bitVec[dependentInd] == ABSENT && bitVec[i]==PRESENT){
                    //throw new OperatorFailedException("");
                    return Double.NEGATIVE_INFINITY;
                }
            }

        }
        parameter.setParameterValue(index, newVal);

        return logq;
    }

    public final String getPerformanceSuggestion() {
        return "no suggestions available";
    }

    public static dr.xml.XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return "msatModelSwitchOperator";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);
            Parameter modelChoose = (Parameter) xo.getElementFirstChild(MODEL_CHOOSE);
            Parameter dependencies = (Parameter)xo.getElementFirstChild(DEPENDENCIES);
            int[] variableIndices;
            if(xo.hasChildNamed(VARIABLE_INDICES)){

                double[] temp = ((Parameter)xo.getElementFirstChild(VARIABLE_INDICES)).getParameterValues();
                variableIndices = new int[temp.length];
                for(int i = 0; i < temp.length;i++){
                    variableIndices[i] = (int)temp[i];
                }

            }else{
                variableIndices = new int[]{0, 1, 2, 3, 4, 5};
            }

            return new MicrosatelliteBitFlipOperator(modelChoose, dependencies, weight, variableIndices);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a microsatellite averaging operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(MODEL_CHOOSE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(DEPENDENCIES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(VARIABLE_INDICES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true)

        };

    };

}
