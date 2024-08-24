/*
 * SequenceErrorModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.tipstatesmodel;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tipstatesmodel.TimeVaryingFrequenciesModel;
import dr.evomodel.treedatalikelihood.TipStateAccessor;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A Suchard
 */
public class TimeVaryingFrequenciesModelParser extends AbstractXMLObjectParser {

    public static final String FREQUENCIES_MODEL = "sequenceErrorModel";
    public static final String CUT_OFF = "cutOff";

    public String getParserName() {
        return FREQUENCIES_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<TipStateAccessor> accessors = xo.getAllChildren(TipStateAccessor.class);

        List<TimeVaryingFrequenciesModel.Epoch> epochs = new ArrayList<>();

        for (XMLObject cxo : xo.getAllChildren(CUT_OFF)) {
            double cutOff = cxo.getDoubleAttribute(ParameterParser.VALUE);
            Parameter freq = (Parameter) cxo.getChild(Parameter.class);

            epochs.add(new TimeVaryingFrequenciesModel.Epoch(cutOff, freq));
        }

        Parameter last = (Parameter) xo.getChild(Parameter.class);
        epochs.add(new TimeVaryingFrequenciesModel.Epoch(Double.POSITIVE_INFINITY, last));

        Taxon taxon = (Taxon) xo.getChild(Taxon.class);
        Tree tree = (Tree) xo.getChild(Tree.class);


        return new TimeVaryingFrequenciesModel(xo.getId(), accessors, epochs, taxon, tree);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a model that allows for post-mortem DNA damage.";
    }

    public Class getReturnType() {
        return TimeVaryingFrequenciesModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CUT_OFF, new XMLSyntaxRule[] {
                    AttributeRule.newDoubleRule(ParameterParser.VALUE),
                    new ElementRule(Parameter.class),
            }, 0, Integer.MAX_VALUE),
            new ElementRule(Parameter.class),
            // TODO Add data-type to check dimensions
    };
}
