/*
 * SpatialKernel.java
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

package dr.evomodel.epidemiology.casetocase;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.IntegrableUnivariateFunction;
import dr.math.RiemannApproximation;
import dr.xml.*;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

/**
 * An abstract spatial kernel function class (and some examples) with a list of parameters.
 *
 * @author Matthew Hall
 */

public abstract class SpatialKernel extends AbstractModel implements IntegrableUnivariateFunction {

    private ArrayList<Parameter> params;

    public SpatialKernel(String name){
        super(name);
    }

    public ArrayList<Parameter> getParams(){
        return params;
    }

    public Parameter getParam(int index){
        return params.get(index);
    }

    public void setParam(int index, Parameter value){
        params.set(index, value);
    }

    public static final String SPATIAL_KERNEL_FUNCTION = "spatialKernelFunction";
    public static final String PARAMETERS = "parameters";
    public static final String KERNEL_TYPE = "type";
    public static final String INTEGRATOR_STEPS = "integratorSteps";
    public static final String ALPHA = "kernel.alpha";
    public static final String R_0 = "kernel.r0";

    public enum Type{
        EXPONENTIAL ("exponential", Exponential.class),
        POWER_LAW ("powerLaw", PowerLaw.class),
        GAUSSIAN ("gaussian", Gaussian.class),
        LOGISTIC ("logistic", Logistic.class);

        private final String xmlName;
        private final Class kernelClass;

        String getXmlName(){
            return xmlName;
        }

        SpatialKernel makeKernelFunction(ArrayList<Parameter> parameters) throws IllegalAccessException,
                InstantiationException, InvocationTargetException, NoSuchMethodException {
//            Constructor[] construct = kernelClass.getConstructors();

            Constructor constructor = kernelClass.getConstructor(SpatialKernel.class, String.class, ArrayList.class);

            return (SpatialKernel)constructor.newInstance(null, xmlName, parameters);
        }

        Type(String xmlName, Class kernelClass){
            this.xmlName = xmlName;
            this.kernelClass = kernelClass;
        }
    }

    public static double EuclideanDistance(double[] point1, double[] point2){
        return Math.sqrt(Math.pow(point1[0]-point2[0],2) + Math.pow(point1[1]-point2[1], 2));
    }

    public double getUpperBound(){
        return Double.POSITIVE_INFINITY;
    }

    public double getLowerBound(){
        return 0;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index){
        fireModelChanged();
    }

    protected void storeState(){
    }

    protected void restoreState(){
    }

    public double value(double distance){
        return evaluate(distance);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type){
        fireModelChanged();
    }

    protected void acceptState(){
        // nothing to do?
    }

    // no need to do this unless there is one...

    public void configureIntegrator(int sampleSize){}

    public abstract SpatialKernel newInstance(ArrayList<Parameter> params) throws InstantiationException;

    public class Exponential extends SpatialKernel  {

        private Parameter alpha;

        public SpatialKernel newInstance(ArrayList<Parameter> params) throws InstantiationException {
            return new Exponential(Type.EXPONENTIAL.getXmlName(), params);
        }

        public Exponential(String name, ArrayList<Parameter> params) throws InstantiationException {
            super(name);

            if(params.size()!=1){
                throw new InstantiationException("Wrong number of parameters for this spatial kernal function");
            }

            if(!params.get(0).getId().equals(ALPHA)){
                throw new InstantiationException("No parameter named alpha");
            }

            this.alpha = params.get(0);
            addVariable(alpha);
        }

        public Exponential() throws InstantiationException {
            this(Type.EXPONENTIAL.getXmlName(), null);
        }



        public double evaluate(double argument) {
            return evaluate(argument, alpha.getParameterValue(0));
        }

        public double evaluate(double argument, double alpha){
            return Math.exp(-argument * alpha);
        }

        public double evaluateIntegral(double a, double b){
            double aValue = alpha.getParameterValue(0);
            return -(1/aValue)*Math.exp(-aValue*b) + (1/aValue)*Math.exp(-aValue*a);
        }

    }

    public class PowerLaw extends SpatialKernel {

        private Parameter alpha;

        public SpatialKernel newInstance(ArrayList<Parameter> params) throws InstantiationException{
            return new PowerLaw(Type.POWER_LAW.getXmlName(), params);
        }

        public PowerLaw(String name, ArrayList<Parameter> params) throws InstantiationException {
            super(name);
            if(params.size()!=1){
                throw new InstantiationException("Wrong number of parameters for this spatial kernal function");
            }

            if(!params.get(0).getId().equals(ALPHA)){
                throw new InstantiationException("No parameter named alpha");
            }
            this.alpha = params.get(0);
            addVariable(alpha);
        }

        public PowerLaw() throws InstantiationException {
            this(Type.POWER_LAW.getXmlName(), null);
        }



        public double evaluate(double argument) {
            return evaluate(argument, alpha.getParameterValue(0));
        }

        public double evaluate(double argument, double alpha){
            return Math.pow(argument, -alpha);
        }

        public double evaluateIntegral(double a, double b){
            double aValue = params.get(0).getParameterValue(0);
            return -aValue*Math.pow(b, -aValue-1) + -aValue*Math.pow(a, -aValue-1);
        }

    }

    public class Gaussian extends SpatialKernel {

        private Parameter alpha;

        RiemannApproximation integrator;

        public SpatialKernel newInstance(ArrayList<Parameter> params) throws InstantiationException{
            return new Gaussian(Type.GAUSSIAN.getXmlName(), params);
        }

        public Gaussian(String name, ArrayList<Parameter> params) throws InstantiationException {
            this(name, params, 25);
        }

        public Gaussian() throws InstantiationException {
            this(Type.GAUSSIAN.getXmlName(), null);
        }

        public void configureIntegrator(int sampleSize){
            integrator = new RiemannApproximation(sampleSize);
        }

        public Gaussian(String name, ArrayList<Parameter> params, int steps) throws InstantiationException {
            super(name);

            if(params.size()!=1){
                throw new InstantiationException("Wrong number of parameters for this spatial kernal function");
            }

            if(!params.get(0).getId().equals(ALPHA)){
                throw new InstantiationException("No parameter named alpha");
            }

            this.alpha = params.get(0);
            addVariable(alpha);


            integrator = new RiemannApproximation(steps);
        }

        public double evaluate(double argument) {
            return evaluate(argument, alpha.getParameterValue(0));
        }

        public double evaluate(double argument, double alpha){
            return Math.exp(-Math.pow(argument, 2) * alpha);
        }

        public double evaluateIntegral(double a, double b){
            return integrator.integrate(this, a, b);
        }

    }

    public class Logistic extends SpatialKernel {

        private Parameter alpha;
        private Parameter r_0;

        RiemannApproximation integrator;

        public SpatialKernel newInstance(ArrayList<Parameter> params) throws InstantiationException{

            return new Logistic(Type.LOGISTIC.getXmlName(), params);
        }

        public Logistic() throws InstantiationException {
            this(Type.GAUSSIAN.getXmlName(), null);
        }

        public Logistic(String name, ArrayList<Parameter> params) throws InstantiationException {
            this(name, params, 25);
        }


        public void configureIntegrator(int sampleSize){
            integrator = new RiemannApproximation(sampleSize);
        }

        public Logistic(String name, ArrayList<Parameter> params, int steps) throws InstantiationException {
            super(name);

            if(params.size()!=2){
                throw new InstantiationException("Wrong number of parameters for this spatial kernal function");
            }

            boolean hasAlpha = false;
            boolean hasR0 = false;

            for(Parameter parameter : params){
                if(parameter.getId().equals(ALPHA)){
                    hasAlpha = true;
                    this.alpha = parameter;
                }
                if( parameter.getId().equals(R_0)){
                    hasR0 = true;
                    this.r_0 = parameter;
                }
            }

            if(!hasAlpha || !hasR0){
                throw new InstantiationException("Kernel function does not have the required parameters");
            }

            addVariable(alpha);
            addVariable(r_0);

            integrator = new RiemannApproximation(steps);
        }

        public double evaluate(double argument) {
            return evaluate(argument, alpha.getParameterValue(0), r_0.getParameterValue(0));
        }

        public double evaluate(double argument, double alpha, double r_0){
            return 1/(1+Math.pow((argument/r_0), alpha));
        }

        public double evaluateIntegral(double a, double b){
            return integrator.integrate(this, a, b);
        }

    }





    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName(){
            return SPATIAL_KERNEL_FUNCTION;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            try{
                String type = (String)xo.getAttribute("type");
                XMLObject parameters = xo.getChild(PARAMETERS);
                ArrayList<Parameter> params = new ArrayList<Parameter>();
                for(int i=0; i<parameters.getChildCount(); i++){
                    Parameter cxo = (Parameter)parameters.getChild(i);
                    params.add(cxo);
                }

                SpatialKernel kernelFunction = null;
                for(Type value : Type.values()){
                    if(value.getXmlName().equals(type)){
                        kernelFunction = value.makeKernelFunction(params);
                    }
                }
                if(kernelFunction==null){
                    throw new XMLParseException("Unknown spatial kernel type");
                }

                if(xo.hasAttribute(INTEGRATOR_STEPS)){
                    kernelFunction.configureIntegrator(Integer.parseInt((String)xo.getAttribute(INTEGRATOR_STEPS)));
                }
                return kernelFunction;
            } catch(Exception e){
                throw new XMLParseException("Failed to initiate spatial kernel ("+e.toString()+")");
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }


        public final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(PARAMETERS, new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
                        }, 1, 1),
                AttributeRule.newStringRule(KERNEL_TYPE),
                AttributeRule.newIntegerRule(INTEGRATOR_STEPS, true)
        };

        public String getParserDescription() {
            return "This element represents a spatial kernel function with a single parameter.";
        }

        public Class getReturnType() {
            return SpatialKernel.class;
        }
    };



}
