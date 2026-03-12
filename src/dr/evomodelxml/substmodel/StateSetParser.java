/*
 * LogCtmcRatesMatrixMatrixProductParameterParser.java
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

package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.StronglyLumpableCtmcRates;
import dr.inference.model.Parameter;
import dr.util.Identifiable;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

//import static dr.evomodel.substmodel.StronglyLumpableCtmcRates.StateSet;
import dr.inference.model.StateSet;
import static dr.evoxml.GeneralDataTypeParser.CODE;
import static dr.evoxml.GeneralDataTypeParser.STATE;

/**
 * @author Xinghua Tao
 * @author Marc A. Suchard
 */

public class StateSetParser extends AbstractXMLObjectParser {

    private static final String STATE_SET = "stateSet";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String id = xo.getId();
        DataType dataType = (DataType) xo.getChild(DataType.class);

        List<Integer> states = new ArrayList<>();
        for (XMLObject cxo : xo.getAllChildren(STATE)) {
            String string = cxo.getStringAttribute(CODE);
            int state = dataType.getState(string);
            if (state < 0) {
                throw new XMLParseException("Unknown code '" + string + "'");
            }
            states.add(state);
        }

        return new StateSet(id, states, dataType);
    }

//    private StateSet parseStateSet(XMLObject xo, DataType dataType) throws XMLParseException {
//
//        }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ContentRule("<state code=\"X\"/>"),
            new ElementRule(Identifiable.class, 1, Integer.MAX_VALUE),
    };

    public String getParserDescription() {
        return "";
    } // TODO

    public Class getReturnType() {
        return StateSet.class;
    }

    public String getParserName() {
        return STATE_SET;
    }
}
