/*
 * MicrosatellitePatternParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evoxml;

import dr.xml.*;
import dr.evolution.util.Taxa;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.alignment.Patterns;
import java.util.logging.Logger;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser that returns a list of microsatellite patterns
 *
 */
public class MicrosatellitePatternParser extends AbstractXMLObjectParser {

    public static final String MICROSATPATTERN = "microsatellitePattern";
    public static final String MICROSAT_SEQ = "microsatSeq";
    public static final String PRINT_DETAILS = "printDetails";
    public static final String PRINT_MSAT_PATTERN_CONTENT = "printMsatPatContent";
    public static final String ID ="id";


    public static final int COUNT_INCREMENT = 100;

    public String getParserName() {
        return MICROSATPATTERN;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Taxa taxonList = (Taxa)xo.getChild(Taxa.class);

        Microsatellite microsatellite = (Microsatellite)xo.getChild(Microsatellite.class);

        String[] strLengths = ((String) xo.getElementFirstChild(MICROSAT_SEQ)).split(",");
        int[] pattern = new int[strLengths.length];
        try {
            for (int i = 0; i < strLengths.length; i++) {
                pattern[i] = microsatellite.getState(strLengths[i]);
            }
        } catch (NumberFormatException nfe) {
            throw new XMLParseException("Unable to parse microsatellite data: " + nfe.getMessage());
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("Unable to parse microsatellite data: " + iae.getMessage());
        }

        Patterns microsatPat = new Patterns(microsatellite, taxonList);
        microsatPat.addPattern(pattern);
        microsatPat.setId((String)xo.getAttribute(ID));

        if(xo.getAttribute(PRINT_DETAILS,true)){
            printDetails(microsatPat);
        }

        if(xo.getAttribute(PRINT_MSAT_PATTERN_CONTENT,true)){
            printMicrosatContent(microsatPat);
        }

        return microsatPat;

    }

    public static void printDetails(Patterns microsatPat){
        Logger.getLogger("dr.evoxml").info(
                "    Locus name: "+microsatPat.getId()+
                        "\n    Number of Taxa: "+microsatPat.getPattern(0).length+
                        "\n    min: "+((Microsatellite)microsatPat.getDataType()).getMin()+" "+
                        "max: "+((Microsatellite)microsatPat.getDataType()).getMax()+
                        "\n    state count: "+microsatPat.getDataType().getStateCount()+"\n");
    }

    public static void printMicrosatContent(Patterns microsatPat){
        Logger.getLogger("dr.evoxml").info(
                "    Locus name: "+ microsatPat.getId());
        int[] pat = microsatPat.getPattern(0);
        for(int i = 0; i < pat.length; i++){
            Logger.getLogger("dr.evoxml").info("    Taxon: "+microsatPat.getTaxon(i)+" "+"state: "+pat[i]);
        }
        Logger.getLogger("dr.evoxml").info("\n");
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Taxa.class),
            new ElementRule(Microsatellite.class),
            new ElementRule(MICROSAT_SEQ,new XMLSyntaxRule[]{
                    new ElementRule(String.class,
                            "A string of numbers representing the microsatellite lengths for a locus",
                            "1,2,3,4,5,67,100")},false),
            new StringAttributeRule(ID, "the name of the locus"),
            AttributeRule.newBooleanRule(PRINT_DETAILS, true),
            AttributeRule.newBooleanRule(PRINT_MSAT_PATTERN_CONTENT, true)

    };

    public String getParserDescription() {
        return "This element represents a microsatellite pattern.";
    }

    public Class getReturnType() {
        return Patterns.class;
    }


}