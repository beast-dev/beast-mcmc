package dr.evomodel.continuous;

import dr.xml.*;

/*
 * @author Marc Suchard
 */
public class BranchDirectionAttributeProvider extends BivariateTraitBranchAttributeProvider {

    public static final String DIRECTION_PROVIDER = "branchDirections";
    public static String TRAIT_EXTENSION = ".angle";

    public BranchDirectionAttributeProvider(SampledMultivariateTraitLikelihood traitLikelihood) {
        super(traitLikelihood);
    }
     
    protected String extensionName() {
        return TRAIT_EXTENSION;
    }

    protected double convert(double latValue, double longValue, double timeValue) {
        double angle = Math.atan2(latValue,longValue);
             if (angle < 0)
                 angle = 2*Math.PI + angle;
        return angle;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             SampledMultivariateTraitLikelihood traitLikelihood = (SampledMultivariateTraitLikelihood)
                     xo.getChild(SampledMultivariateTraitLikelihood.class);

             return new BranchDirectionAttributeProvider(traitLikelihood);
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return new XMLSyntaxRule[] {
                 new ElementRule(SampledMultivariateTraitLikelihood.class),
             };
         }

         public String getParserDescription() {
             return null;
         }

         public Class getReturnType() {
             return BivariateTraitBranchAttributeProvider.class;
         }

         public String getParserName() {
             return DIRECTION_PROVIDER;
         }
     };

//     public static void main(String[] arg) {
//        double[] vector = {-10,1};
//
//        double angle = convert(vector[0],vector[1]);
//
//        System.err.println("vec:   "+new Vector(vector));
//        System.err.println("angle: "+angle);
//    }

}
