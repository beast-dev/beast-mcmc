package dr.evomodel.epidemiology.casetocase;

import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * An abstract spatial kernel function class (and some examples) with a single parameter a.
 *
 * @author Matthew Hall
 */

public abstract class SpatialKernelFunction {

    private Parameter a;

    public SpatialKernelFunction(Parameter a){
        this.a = a;
    }

    public Parameter geta(){
        return a;
    }

    public void seta(Parameter value){
        a = value;
    }

    public abstract double value(double[] point1, double[] point2);

    public static final String SPATIAL_KERNEL_FUNCTION = "spatialKernelFunction";
    public static final String A = "a";
    public static final String KERNEL_TYPE = "kernelType";

    public enum Type{
        EXPONENTIAL ("exponential", Exponential.class),
        POWER_LAW ("powerLaw", PowerLaw.class),
        GAUSSIAN ("gaussian", Gaussian.class);

        private final String xmlName;
        private final Class kernelClass;

        String getXmlName(){
            return xmlName;
        }

        SpatialKernelFunction makeKernelFunction() throws InstantiationException, IllegalAccessException{
            SpatialKernelFunction out = (SpatialKernelFunction)kernelClass.newInstance();
            return out;
        }

        Type(String xmlName, Class kernelClass){
            this.xmlName = xmlName;
            this.kernelClass = kernelClass;
        }
    }

    public static double EuclideanDistance(double[] point1, double[] point2){
        return Math.sqrt(Math.pow(point1[0]-point2[0],2) + Math.pow(point1[1]-point2[1], 2));
    }

    public class Exponential extends SpatialKernelFunction {

        private Parameter a;

        public Exponential(Parameter a){
            super(a);
        }

        public double value(double[] point1, double[] point2){
            return Math.exp(-EuclideanDistance(point1, point2)*a.getParameterValue(0));
        }

    }

    public class PowerLaw extends SpatialKernelFunction {

        Parameter a;

        public PowerLaw(Parameter a){
            super(a);
        }

        public double value(double[] point1, double[] point2){
            return Math.pow(EuclideanDistance(point1, point2),-a.getParameterValue(0));
        }


    }

    public class Gaussian extends SpatialKernelFunction {

        private Parameter a;

        public Gaussian(Parameter a){
            super(a);
        }

        public double value(double[] point1, double[] point2){
            return Math.exp(-Math.pow(EuclideanDistance(point1, point2),2)*a.getParameterValue(0));
        }

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName(){
            return SPATIAL_KERNEL_FUNCTION;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            try{
                String type = (String)xo.getAttribute("type");
                SpatialKernelFunction kernelFunction = null;
                for(Type value : Type.values()){
                    if(value.getXmlName().equals(type)){
                        kernelFunction = value.makeKernelFunction();
                    }
                }
                if(kernelFunction==null){
                    throw new XMLParseException("Unknown spatial kernel type");
                }
                Parameter a = new Parameter.Default((Double) xo.getAttribute(A));
                kernelFunction.seta(a);
                return kernelFunction;
            } catch(InstantiationException e){
                throw new XMLParseException("Failed to initiate spatial kernel");
            } catch(IllegalAccessException e){
                throw new XMLParseException("Failed to initiate spatial kernel");
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }


        public final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(A),
                    AttributeRule.newStringRule(KERNEL_TYPE)
            };
        }

        public String getParserDescription() {
            return "This element represents a spatial kernel function with a single parameter.";
        }

        public Class getReturnType() {
            return SpatialKernelFunction.class;
        }
    };



}
