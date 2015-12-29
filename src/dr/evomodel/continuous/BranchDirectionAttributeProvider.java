/*
 * BranchDirectionAttributeProvider.java
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
