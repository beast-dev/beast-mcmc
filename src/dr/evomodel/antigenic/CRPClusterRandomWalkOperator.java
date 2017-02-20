/*
 * CRPClusterRandomWalkOperator.java
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

package dr.evomodel.antigenic;



import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.logging.Logger;


import dr.inference.model.Parameter;


/**
 * Created with IntelliJ IDEA.
 * User: Gabi
 * Date: 2/2/15
 * Time: 11:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class CRPClusterRandomWalkOperator extends AbstractCoercableOperator {



    public static final String WINDOW_SIZE = "windowSize";
    public static final String VIRUS_LOCATIONS = "virusLocations";
    public static final String CLUSTER_ASSIGNMENTS = "clusterAssignments";
    public static final String CRP_CLUSTER_RW_OPERATOR= "CRPClusterRandomWalkOperator";



    private final MatrixParameter virusLocations;
    private final Parameter assignments;
    private double windowSize;

    public CRPClusterRandomWalkOperator( Parameter assignments, MatrixParameter virusLocations,double windowSize,   double weight){
        super(CoercionMode.COERCION_ON);

        this.virusLocations=virusLocations;
        this.windowSize = windowSize;
        this.assignments=assignments;
        setWeight(weight);


    }


    public double doOperation() {

       int max = findMax(assignments);




       int group = MathUtils.nextInt(max+1);


        if (windowSize>50){
            windowSize=50;
        }

        double draw0 = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;
        double draw1 = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;


       for (int i = 0; i<assignments.getDimension(); i++){
           if ((int)assignments.getParameterValue(i)==group){
               double newValue0 = virusLocations.getParameter(i).getParameterValue(0) + draw0;
               double newValue1 = virusLocations.getParameter(i).getParameterValue(1) + draw1;


               virusLocations.setParameterValue(2*i,newValue0);
               virusLocations.setParameterValue(2*i+1,newValue1);

               virusLocations.fireParameterChangedEvent(2*i,Variable.ChangeType.VALUE_CHANGED);
               virusLocations.fireParameterChangedEvent(2*i+1,Variable.ChangeType.VALUE_CHANGED);

               //   virusLocations.getParameter(i).setParameterValue(0,newValue0);
           }
       }


    //   printInformation(virusLocations, "depois");



        return 0.0;
    }





    private int findMax(Parameter assignments){
        int max=0;
        for (int i=0; i<assignments.getDimension(); i++ ){
           if (assignments.getParameterValue(i)>max){
               max = (int) assignments.getParameterValue(i);
           }
        }
    return max;
    }




    private void printInformation(double[] par, int length) {
        StringBuffer sb = new StringBuffer("\n \n double vector \n");
        for(int j=0; j<length; j++){
            sb.append(par[j] + "\t");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString()); };





    private void printInformation(Parameter par, String lala) {
        StringBuffer sb = new StringBuffer("\n ");
        sb.append(lala);
        sb.append("\t\t") ;

        for (int i=0; i<par.getDimension(); i++){
        sb.append(par.getParameterValue(i)+ "\t");                           }

        Logger.getLogger("dr.evomodel").info(sb.toString()); };





    private void printInformation(int par, String lala) {
        StringBuffer sb = new StringBuffer("\n");

        sb.append(lala);
        sb.append("\t\t");
        sb.append(par);

        Logger.getLogger("dr.evomodel").info(sb.toString()); };





    private void printInformation(double par, String lala) {
        StringBuffer sb = new StringBuffer("\n");

        sb.append(lala);
        sb.append("\t\t");
        sb.append(par);

        Logger.getLogger("dr.evomodel").info(sb.toString()); };




    public String getPerformanceSuggestion() {
        return null;
    }



    public String getOperatorName() {
        return CRP_CLUSTER_RW_OPERATOR;
    }



    public double getCoercableParameter() {
        return Math.log(windowSize);
    }

    public void setCoercableParameter(double value) {
        windowSize = Math.exp(value);
    }

    public double getRawParameter() {
        return windowSize;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }




    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return CRP_CLUSTER_RW_OPERATOR;
        }



        public Object parseXMLObject(XMLObject xo) throws XMLParseException {




            double weight = xo.getDoubleAttribute(WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);


            MatrixParameter virusLocations = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
            Parameter assignments = (Parameter) xo.getElementFirstChild(CLUSTER_ASSIGNMENTS);
            return new CRPClusterRandomWalkOperator(assignments,virusLocations,windowSize,weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an operator moves simultaneously locations for all viruses in same cluster in ddCRP and BMDS.";
        }

        public Class getReturnType() {
            return CRPClusterRandomWalkOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class),
                new ElementRule(CLUSTER_ASSIGNMENTS, Parameter.class),

                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newDoubleRule(WINDOW_SIZE),
        };
    };
}






