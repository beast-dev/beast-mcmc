/*
 * GMRFSkyrideLikelihoodParser.java
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

package dr.evomodelxml.coalescent;


import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.GMRFSkyrideGradient;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.hmc.GMRFGradient;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightTransform;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class GMRFSkyrideGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "gmrfSkyrideGradient";
    private static final String WRT_PARAMETER = "wrtParameter";

    private static final String COALESCENT_INTERVAL = "coalescentInterval";
    private static final String NODE_HEIGHT = "nodeHeight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        GMRFSkyrideLikelihood skyrideLikelihood = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);

        String wrtParameterCase = (String) xo.getAttribute(WRT_PARAMETER);

        GMRFGradient.WrtParameter type = GMRFGradient.WrtParameter.factory(wrtParameterCase);
        if (type != null) {
            return new GMRFGradient((GMRFMultilocusSkyrideLikelihood) skyrideLikelihood, type);
        }

        // Old behaviour

        GMRFSkyrideGradient.WrtParameter wrtParameter = setupWrtParameter(wrtParameterCase);

        NodeHeightTransform nodeHeightTransform = (NodeHeightTransform) xo.getChild(NodeHeightTransform.class);

        return new GMRFSkyrideGradient(skyrideLikelihood, wrtParameter, parameter, nodeHeightTransform);
    }

    private GMRFSkyrideGradient.WrtParameter setupWrtParameter(String wrtParameterCase) {
        if (wrtParameterCase.equalsIgnoreCase(COALESCENT_INTERVAL)) {
            return GMRFSkyrideGradient.WrtParameter.COALESCENT_INTERVAL;
        } else if (wrtParameterCase.equalsIgnoreCase(NODE_HEIGHT)) {
            return GMRFSkyrideGradient.WrtParameter.NODE_HEIGHTS;
        }
        else {
            throw new RuntimeException("Not yet implemented!");
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(WRT_PARAMETER),
            new ElementRule(Parameter.class, true),
            new ElementRule(GMRFSkyrideLikelihood.class),
            new ElementRule(NodeHeightTransform.class, true),
    };

    @Override
    public Class getReturnType() {
        return GMRFSkyrideGradient.class;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
