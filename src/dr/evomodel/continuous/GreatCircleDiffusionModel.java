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
//    public static final String BIAS = "mu";
//    public static final String PRECISION_TREE_ATTRIBUTE = "precision";

    public static final double LOG2PI = Math.log(2*Math.PI);

    public GreatCircleDiffusionModel(Parameter precision) {
        super();
        this.precision = precision;
        addParameter(precision);
    }


    public double getLogLikelihood(double[] start, double[] stop, double time) {

        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(start[0], start[1]);
        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(stop[0],  stop[1] );

        double distance = coord1.distance(coord2);

        double prec = precision.getParameterValue(0);
        return 0.5 * (Math.log(prec) - LOG2PI - distance * distance * prec);
    }

    protected void calculatePrecisionInfo() {
    }    


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

         public String getParserName() {
             return DIFFUSION_PROCESS;
         }

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//             XMLObject cxo = (XMLObject) xo.getChild(DIFFUSION_CONSTANT);
             Parameter diffusionParam = (Parameter) xo.getChild(Parameter.class);

             return new GreatCircleDiffusionModel(diffusionParam);
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
//                 new ElementRule(DIFFUSION_CONSTANT,
//                         new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
         };

         public Class getReturnType() {
             return MultivariateDiffusionModel.class;
         }
     };

    private Parameter precision;


}
