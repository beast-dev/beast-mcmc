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

/**
 * An abstract spatial kernel function class (and some examples) with single parameter aParam.
 *
 * @author Matthew Hall
 */

public abstract class SpatialKernel extends AbstractModel implements IntegrableUnivariateFunction {

    private Parameter aParam;

    public SpatialKernel(String name, Parameter a){
        super(name);
        this.aParam = a;
        this.addVariable(a);
    }

    public Parameter geta(){
        return aParam;
    }

    public void seta(Parameter value){
        aParam = value;
    }

    public static final String SPATIAL_KERNEL_FUNCTION = "spatialKernelFunction";
    public static final String A = "a";
    public static final String KERNEL_TYPE = "type";
    public static final String INTEGRATOR_STEPS = "integratorSteps";

    public enum Type{
        EXPONENTIAL ("exponential", Exponential.class),
        POWER_LAW ("powerLaw", PowerLaw.class),
        GAUSSIAN ("gaussian", Gaussian.class);

        private final String xmlName;
        private final Class kernelClass;

        String getXmlName(){
            return xmlName;
        }

        SpatialKernel makeKernelFunction(Parameter a) throws IllegalAccessException, InstantiationException,
                InvocationTargetException {
            Constructor[] construct = kernelClass.getConstructors();
            return (SpatialKernel)construct[0].newInstance(null, xmlName, a);
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

    public double value(double distance, double alpha){
        return evaluate(distance, alpha);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type){
        fireModelChanged();
    }

    protected void acceptState(){
        // nothing to do?
    }

    public double evaluate(double argument){
        return evaluate(argument, aParam.getParameterValue(0));
    }

    // no need to do this unless there is one...

    public void configureIntegrator(int sampleSize){}

    public abstract double evaluate(double argument, double alpha);

    public abstract SpatialKernel newInstance(Parameter a);

    public class Exponential extends SpatialKernel {

        public SpatialKernel newInstance(Parameter a){
            return new Exponential(Type.EXPONENTIAL.getXmlName(), a);
        }

        public Exponential(){
            this(Type.EXPONENTIAL.getXmlName(), null);
        }

        public Exponential(String name, Parameter a){
            super(name, a);
        }

        public double evaluate(double argument, double alpha){
            return Math.exp(-argument * alpha);
        }

        public double evaluateIntegral(double a, double b){
            double aValue = aParam.getParameterValue(0);
            return -(1/aValue)*Math.exp(-aValue*b) + (1/aValue)*Math.exp(-aValue*a);
        }

    }

    public class PowerLaw extends SpatialKernel {

        public SpatialKernel newInstance(Parameter a){
            return new PowerLaw(Type.POWER_LAW.getXmlName(), a);
        }

        public PowerLaw(){
            this(Type.POWER_LAW.getXmlName(), null);
        }

        public PowerLaw(String name, Parameter a){
            super(name, a);
        }

        public double value(double[] point1, double[] point2){
            return evaluate(EuclideanDistance(point1, point2));
        }

        public double evaluate(double argument, double alpha){
            return Math.pow(argument, -alpha);
        }

        public double evaluateIntegral(double a, double b){
            double aValue = aParam.getParameterValue(0);
            return -aValue*Math.pow(b, -aValue-1) + -aValue*Math.pow(a, -aValue-1);
        }

    }

    public class Gaussian extends SpatialKernel {

        RiemannApproximation integrator;

        public SpatialKernel newInstance(Parameter a){
            return new Gaussian(Type.GAUSSIAN.getXmlName(), a);
        }

        public Gaussian(){
            this(Type.GAUSSIAN.getXmlName(), null);
        }

        public Gaussian(String name, Parameter a){
            this(name, a, 25);
        }

        public void configureIntegrator(int sampleSize){
            integrator = new RiemannApproximation(sampleSize);
        }

        public Gaussian(String name, Parameter a, int steps){
            super(name, a);
            integrator = new RiemannApproximation(steps);
        }

        public double evaluate(double argument, double alpha){
            return Math.exp(-Math.pow(argument, 2) * alpha);
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
                Parameter a = (Parameter)xo.getElementFirstChild(A);
                SpatialKernel kernelFunction = null;
                for(Type value : Type.values()){
                    if(value.getXmlName().equals(type)){
                        kernelFunction = value.makeKernelFunction(a);
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


        public final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    new ElementRule(A, Parameter.class, "The single parameter of this kernel"),
                    AttributeRule.newStringRule(KERNEL_TYPE),
                    AttributeRule.newIntegerRule(INTEGRATOR_STEPS, true)
            };
        }

        public String getParserDescription() {
            return "This element represents a spatial kernel function with a single parameter.";
        }

        public Class getReturnType() {
            return SpatialKernel.class;
        }
    };



}
