/*
 * GeneralF81ModelParser.java
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

import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.GeneralF81Model;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Parses a GeneralF81Model
 */
public class GeneralF81ModelParser extends AbstractXMLObjectParser {

    public static final String GENERAL_F81_MODEL = "generalF81Model";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return GENERAL_F81_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        FrequencyModel freqModel = null;
        if (xo.hasChildNamed(FREQUENCIES)) {
            XMLObject cxo = xo.getChild(FREQUENCIES);
            freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        }
        Logger.getLogger("dr.evomodel").info("  General F81 Model from frequencyModel '"+ freqModel.getId() + "' (stateCount=" + freqModel.getFrequencyCount() + ")");

        return new GeneralF81Model(freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general F81 model for an arbitrary number of states.";
    }

    public Class getReturnType() {
        return GeneralF81ModelParser.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES, FrequencyModel.class)
    };
}
