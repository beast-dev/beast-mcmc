/*
 * BranchMagnitudeAttributeProvider.java
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

import dr.xml.*;

/*
 * @author Marc Suchard
 */
public class BranchMagnitudeAttributeProvider extends BivariateTraitBranchAttributeProvider {

    public static final String MAGNITUDE_PROVIDER = "branchMagnitudes";
    public static final String SCALE = "scaleByLength";
    public static String DISTANCE_EXTENSION = ".distance";
    public static final String RATE_EXTENSION = ".rate";

    public BranchMagnitudeAttributeProvider(AbstractMultivariateTraitLikelihood traitLikelihood, boolean scale) {
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

             AbstractMultivariateTraitLikelihood traitLikelihood = (AbstractMultivariateTraitLikelihood)
                     xo.getChild(SampledMultivariateTraitLikelihood.class);

             boolean scale = xo.getAttribute(SCALE,false);

             return new BranchMagnitudeAttributeProvider(traitLikelihood, scale);
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return new XMLSyntaxRule[] {
                 new ElementRule(AbstractMultivariateTraitLikelihood.class),
                 AttributeRule.newBooleanRule(SCALE,true),
             };
         }

         public String getParserDescription() {
             return null;
         }

         public Class getReturnType() {
             return BranchMagnitudeAttributeProvider.class;
         }

         public String getParserName() {
             return MAGNITUDE_PROVIDER;
         }
     };

    private boolean scale;

}