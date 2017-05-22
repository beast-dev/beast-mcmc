/*
 * SampleStateModelParser.java
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

package dr.oldevomodelxml.sitemodel;

import dr.oldevomodel.sitemodel.SampleStateModel;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Vector;

/**
 */
public class SampleStateModelParser extends AbstractXMLObjectParser {

    public static final String SAMPLE_STATE_MODEL = "sampleStateModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String PROPORTIONS = "proportions";

    public String getParserName() {
        return SAMPLE_STATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MUTATION_RATE);
        Parameter muParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PROPORTIONS);
        Parameter proportionParameter = (Parameter) cxo.getChild(Parameter.class);

        Vector<Object> subModels = new Vector<Object>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof SubstitutionModel) {
                subModels.addElement(xo.getChild(i));
            }
        }

        return new SampleStateModel(muParam, proportionParameter, subModels);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteModel that has a discrete distribution of substitution models over sites, " +
                "designed for sampling of internal states.";
    }

    public Class getReturnType() {
        return SampleStateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MUTATION_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PROPORTIONS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE)
    };
}
