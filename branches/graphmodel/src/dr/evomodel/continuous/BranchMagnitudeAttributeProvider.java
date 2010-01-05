package dr.evomodel.continuous;

import dr.xml.*;
import dr.evolution.tree.BranchAttributeProvider;

/*
 * @author Marc Suchard
 */
public class BranchMagnitudeAttributeProvider extends BivariateTraitBranchAttributeProvider {

    public static final String MAGNITUDE_PROVIDER = "branchMagnitudes";
    public static final String SCALE = "scaleByLength";
    public static String DISTANCE_EXTENSION = ".distance";
    public static final String RATE_EXTENSION = ".rate";

    public BranchMagnitudeAttributeProvider(SampledMultivariateTraitLikelihood traitLikelihood, boolean scale) {
        super(traitLikelihood);
        this.scale = scale;
        label = traitName + extensionName(); 
    }

    protected String extensionName() {
        if (scale)
            return RATE_EXTENSION;
        return DISTANCE_EXTENSION;
    }

    protected double convert(double latValue, double longValue, double timeValue) {
        double result = Math.sqrt(latValue*latValue + longValue*longValue);
        if (scale)
            result /= timeValue;
        return result;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             SampledMultivariateTraitLikelihood traitLikelihood = (SampledMultivariateTraitLikelihood)
                     xo.getChild(SampledMultivariateTraitLikelihood.class);

             boolean scale = xo.getAttribute(SCALE,false);

             return new BranchMagnitudeAttributeProvider(traitLikelihood, scale);
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return new XMLSyntaxRule[] {
                 new ElementRule(SampledMultivariateTraitLikelihood.class),
                 AttributeRule.newBooleanRule(SCALE,true),
             };
         }

         public String getParserDescription() {
             return null;
         }

         public Class getReturnType() {
             return BranchAttributeProvider.class;
         }

         public String getParserName() {
             return MAGNITUDE_PROVIDER;
         }
     };

    private boolean scale;

}