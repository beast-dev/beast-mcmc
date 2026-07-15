/*
 * PredictiveBeagleSequenceSimulatorParser.java
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

package dr.evomodelxml.treedatalikelihood.pps;

import dr.evolution.alignment.PatternList;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.pps.PredictiveBeagleSequenceSimulatorGenerator;
import dr.evomodel.treedatalikelihood.pps.PredictiveDataGenerator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses a predictiveBeagleSequenceSimulator element
 */
public class PredictiveBeagleSequenceSimulatorParser extends AbstractXMLObjectParser {

    public static final String PREDICTIVE_BEAGLE_SEQUENCE_SIMULATOR = "predictiveBeagleSequenceSimulator";

    public String getParserName() {
        return PREDICTIVE_BEAGLE_SEQUENCE_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        try {
            return new PredictiveBeagleSequenceSimulatorGenerator(patterns, treeDataLikelihood);
        } catch (RuntimeException e) {
            throw new XMLParseException("Error constructing " + PREDICTIVE_BEAGLE_SEQUENCE_SIMULATOR +
                    " element '" + xo.getId() + "': " + e.getMessage());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PatternList.class, false),
            new ElementRule(TreeDataLikelihood.class, false)
    };

    public String getParserDescription() {
        return "For use in posteriorPredictiveLogger. Wraps model information and data.";
    }

    public Class getReturnType() {
        return PredictiveDataGenerator.class;
    }
}
