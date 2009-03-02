package dr.evomodel.continuous;

import dr.xml.*;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.evolution.continuous.SphericalPolarCoordinates;

/**
 * @author Marc A. Suchard
 */

public class GreatCircleDiffusionModel extends MultivariateDiffusionModel {

    public static final String DIFFUSION_PROCESS = "greatCircleDiffusionModel";
    public static final String DIFFUSION_CONSTANT = "precision";
    public static final String COEFFICIENT = "diffusionCoefficient";
//    public static final String BIAS = "mu";
//    public static final String PRECISION_TREE_ATTRIBUTE = "precision";

    public static final double LOG2PI = Math.log(2*Math.PI);

    public GreatCircleDiffusionModel(Parameter precision, Parameter coefficient) {
        super();
        this.precision = precision;
        addParameter(precision);
        this.coefficient = coefficient;
        if (coefficient != null)
            addParameter(coefficient);
    }

    public GreatCircleDiffusionModel(Parameter precision) {
        this(precision,null);
    }

    public double getLogLikelihood(double[] start, double[] stop, double time) {

        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(start[0], start[1]);
        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(stop[0],  stop[1] );

        double distance = coord1.distance(coord2);

        double inverseVariance = precision.getParameterValue(0) / time;
        if (coefficient == null)
            return 0.5 * (Math.log(inverseVariance) - LOG2PI - distance * distance * inverseVariance);
        double coef = - coefficient.getParameterValue(0);
        return 0.5 * (coef * Math.log(inverseVariance) - LOG2PI - distance * distance * Math.pow(inverseVariance,coef));
    }

    protected void calculatePrecisionInfo() {
    }    


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

         public String getParserName() {
             return DIFFUSION_PROCESS;
         }

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             Parameter diffusionParam = (Parameter) xo.getChild(Parameter.class);
             Parameter coefficient = null;
             if (xo.hasChildNamed(COEFFICIENT))
                coefficient = (Parameter) ((XMLObject)xo.getChild(COEFFICIENT)).getChild(Parameter.class);

             return new GreatCircleDiffusionModel(diffusionParam, coefficient);
         }

         //************************************************************************
         // AbstractXMLObjectParser implementation
         //************************************************************************

         public String getParserDescription() {
             return "Describes a bivariate diffusion process using great circle distances.";
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return rules;
         }

         private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                 new ElementRule(Parameter.class),
                 new ElementRule(COEFFICIENT,
                         new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
         };

         public Class getReturnType() {
             return MultivariateDiffusionModel.class;
         }
     };

    private Parameter precision;
    private Parameter coefficient;

}
