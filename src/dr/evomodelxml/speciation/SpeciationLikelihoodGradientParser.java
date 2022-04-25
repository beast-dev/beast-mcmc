/*
 * SpeciationLikelihoodGradientParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class SpeciationLikelihoodGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "speciationLikelihoodGradient";
    private static final String WRT_PARAMETER = "wrtParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SpeciationLikelihood likelihood = (SpeciationLikelihood) xo.getChild(SpeciationLikelihood.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        String wrtParamter = (String) xo.getAttribute(WRT_PARAMETER);

        if (! (likelihood.getSpeciationModel() instanceof BirthDeathGernhard08Model)) {
            throw new RuntimeException("Not yet implemented for other cases.");
        }

        SpeciationLikelihoodGradient.WrtParameter type = SpeciationLikelihoodGradient.WrtParameter.factory(wrtParamter);
        return new SpeciationLikelihoodGradient(likelihood, tree, type);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(WRT_PARAMETER),
            new ElementRule(SpeciationLikelihood.class),
            new ElementRule(TreeModel.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SpeciationLikelihoodGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
