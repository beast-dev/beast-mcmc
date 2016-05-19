/*
 * AttributeUncertainPatternsParser.java
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

package dr.evoxml;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleSiteList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evoxml.util.DataTypeUtils;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class UncertainAttributePatternsParser extends AbstractXMLObjectParser {

    public static final String ATTRIBUTE = AttributePatternsParser.ATTRIBUTE;
    public static final String NAME = "uncertainAttributePatterns";
    public static final String LOCATION_TOKEN = "\\s";
    public static final String PROBABILITY_TOKEN = ":";
    public static final String NORMALIZE = "normalize";

    public String getParserName() { return NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String attributeName = xo.getStringAttribute(ATTRIBUTE);
        TaxonList taxa = (TaxonList) xo.getChild(TaxonList.class);
        DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) {
            throw new XMLParseException("dataType expected for attributePatterns element");
        }

        // using a SimpleSiteList rather than Patterns to allow ancestral reconstruction
        UncertainSiteList patterns = new UncertainSiteList(dataType, taxa);

        boolean normalize = xo.getAttribute(NORMALIZE, true);

        if (dataType == null) { // TODO Is this necessary given XMLSyntaxRules?
            throw new XMLParseException("dataType expected for attributePatterns element");
        }

        double[][] uncertainPattern = new double[taxa.getTaxonCount()][];

        // Parse attributes
        boolean attributeFound = false;

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);

            Object value = taxon.getAttribute(attributeName);

            if (value != null) {
                attributeFound = true;
                List<StateProbability> stateProbabilities;

                try {
                    stateProbabilities = parseStates(value.toString(), dataType);
                } catch (StateParseException e) {
                    throw new XMLParseException("State or probability for attribute (" + attributeName + ") in taxon "
                            + taxon.getId()  + " is invalid; state = \"" + e.getState() + "\" and probability =\""
                            + e.getProbability() + "\"");
                }

                uncertainPattern[i] = convertToPartials(stateProbabilities, dataType, normalize);
            } else {
                throw new XMLParseException("State for attribute (" + attributeName + ") in taxon "
                        + taxon.getId() + " is unknown.");
            }
        }

        if (!attributeFound) {
            throw new XMLParseException("The attribute (" + attributeName + ") was missing in all taxa. Check the name of the attribute.");
        }

        patterns.addPattern(uncertainPattern);

        Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating an uncertain attribute model for attribute \""
                + attributeName + "\"");
        Logger.getLogger("dr.evolution").info("\tIf you publish results using this model, please reference:");
        Logger.getLogger("dr.evolution").info("\t" + Citable.Utils.getCitationString(patterns));
        Logger.getLogger("dr.evolution").info("\n");

        return patterns;
    }


    class StateProbability {
        int state;
        double probability;

        public StateProbability(int state, double probability) {
            this.state = state;
            this.probability = probability;
        }

        public int getState() { return state; }

        public double getProbability() { return probability; }
    }

    class StateParseException extends Exception {
        String state;
        String probability;

        public StateParseException(String state, String probability) {
            this.state = state;
            this.probability = probability;
        }

        public String getState() { return state; }

        public String getProbability() { return probability; }
    }

    private List<StateProbability> parseStates(String string, DataType dataType)
            throws StateParseException {

        List<StateProbability> stateProbabilities = new ArrayList<StateProbability>();

        String[] tokens = string.split(LOCATION_TOKEN);
        for (String token : tokens) {
            String[] component = token.split(PROBABILITY_TOKEN);

            int state = dataType.getState(component[0]);

            double probability = 1.0;

            if (component.length > 1) {
                try {
                    probability = Double.valueOf(component[1]);
                } catch (NumberFormatException e) {
                    probability = Double.NaN;
                }
            }

            if (state < 0 || Double.isNaN(probability) || probability <= 0.0 || probability > 1.0) {
                throw new StateParseException(component[0], (component.length == 1 ? "" : component[1]));
            }

            stateProbabilities.add(new StateProbability(state, probability));
        }

        return stateProbabilities;
    }

    private void normalize(double[] vec) {
        double sum = 0.0;
        for (double x : vec) {
            sum += x;
        }
        for (int i = 0; i < vec.length; ++i) {
            vec[i] /= sum;
        }
    }

    private double[] convertToPartials(List<StateProbability> stateProbabilities, DataType dataType,
                                       boolean normalize) {
        double[] partials = new double[dataType.getStateCount()];

        for (StateProbability state : stateProbabilities) {
            partials[state.getState()] = state.getProbability();
        }

        if (normalize) {
            normalize(partials);
        }
        return partials;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new XORRule(
                new StringAttributeRule(
                    DataType.DATA_TYPE,
                    "The data type",
                    DataType.getRegisteredDataTypeNames(), false),
                new ElementRule(DataType.class)
            ),
            AttributeRule.newStringRule(ATTRIBUTE),
            AttributeRule.newBooleanRule(NORMALIZE, true),
            new ElementRule(TaxonList.class, "The taxon set")
    };

    public String getParserDescription() {
        return "A site pattern defined by an attribute in a set of taxa.";
    }

    public Class getReturnType() { return PatternList.class; }

}