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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static dr.evomodel.substmodel.StronglyLumpableCtmcRates.StateSet;
import static dr.evomodel.substmodel.StronglyLumpableCtmcRates.Proportion;
import static dr.evomodel.substmodel.StronglyLumpableCtmcRates.Lump;

import static dr.evoxml.GeneralDataTypeParser.CODE;
import static dr.evoxml.GeneralDataTypeParser.STATE;

/**
 * @author Xinghua Tao
 * @author Marc A. Suchard
 */

public class StronglyLumpableCtmcRatesParser extends AbstractXMLObjectParser {

    private static final String LUMPABLE_MODEL = "stronglyLumpableCtmcRates";
    private static final String LUMP = "lump";
    private static final String RATES = "rates";
    private static final String PROPORTIONS = "proportions";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Lump> lumps = new ArrayList<>();

        DataType dataType = (DataType) xo.getChild(DataType.class);

        for (XMLObject cxo : xo.getAllChildren(LUMP)) {
            lumps.add(parseLump(cxo, dataType));
        }

        if (lumps.size() < 2) {
            throw new XMLParseException("Must specify at least 2 lumps");
        }

        Parameter rates = (Parameter) xo.getChild(Parameter.class);

        String name = xo.hasId() ? xo.getId() : LUMPABLE_MODEL;

        // TODO Check (and set) dimensions
        return new StronglyLumpableCtmcRates(name, lumps, rates, dataType);
    }

    private Lump parseLump(XMLObject xo, DataType dataType) throws XMLParseException {

        String id = xo.getId();

        StateSet set = (StateSet) xo.getChild(StateSet.class);

        Parameter rates = null;
        if (xo.hasChildNamed(RATES)) {
            rates = (Parameter) xo.getElementFirstChild(RATES);
        } else if (set.size() > 1) {
            throw new XMLParseException("Must provide within-lump rate parameter for '" + id + "'");
        }

        List<Proportion> proportions = new ArrayList<>();
        for (XMLObject cxo : xo.getAllChildren(PROPORTIONS)) {
            String string = cxo.getChild(STATE).getStringAttribute(CODE);
            int source = dataType.getState(string);
            if (source < 0) {
                throw new XMLParseException("Unknown code '" + string + "'");
            }
            StateSet destination = (StateSet) cxo.getChild(StateSet.class);
            Parameter parameter = (Parameter) cxo.getChild(Parameter.class);

            int size = destination.size();
            parameter.setDimension(size);
            for (int i = 0; i < size; ++i) {
                parameter.setParameterValue(i, 1.0 / (double) size); // TODO remove magic values
            }

            proportions.add(new Proportion(source, destination, parameter));
        }

        // TODO Add in destination proportions of size == 1 automatically



        return new Lump(id, set.states(), rates, proportions);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(LUMP, new XMLSyntaxRule[]{
                    new ElementRule(Identifiable.class, 0, Integer.MAX_VALUE),
                    new ContentRule("<state code=\"X\"/>"),
                    new ElementRule(RATES, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
                    new ElementRule(PROPORTIONS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            }, 1, Integer.MAX_VALUE),
            new ElementRule(Parameter.class),
    };

    public String getParserDescription() {
        return "A matrix-matrix product of parameters for CTMC rates.";
    }

    public Class getReturnType() {
        return StronglyLumpableCtmcRates.class;
    }

    public String getParserName() {
        return LUMPABLE_MODEL;
    }
}
