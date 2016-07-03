/*
 * NewMicrosatelliteModelParser.java
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

package dr.oldevomodelxml.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.NewMicrosatelliteModel;
import dr.xml.*;


/**
 * @author Chieh-Hsi Wu
 *
 */
public class NewMicrosatelliteModelParser extends AbstractXMLObjectParser{
    public static final String NEW_MSAT_MODEL = "newMsatModel";

    public String getParserName() {
        return NEW_MSAT_MODEL;
    }


    //AbstractXMLObjectParser implementation
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        //get msat data type
        System.out.println("Using watkins' model");
        Microsatellite msat = (Microsatellite)xo.getChild(Microsatellite.class);
        //get FrequencyModel
        FrequencyModel freqModel = null;
        if(xo.hasChildNamed(FrequencyModelParser.FREQUENCIES)){
            freqModel = (FrequencyModel)xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        }


        return new NewMicrosatelliteModel(msat, freqModel);
    }


    public String getParserDescription() {
        return "This element represents an instance of the stepwise mutation model of microsatellite evolution.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Microsatellite.class),
            new ElementRule(FrequencyModel.class,true)
    };

    public Class getReturnType() {
        return NewMicrosatelliteModel.class;
    }
}
