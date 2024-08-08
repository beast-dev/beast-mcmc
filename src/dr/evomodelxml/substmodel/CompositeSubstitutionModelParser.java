/*
 * GeneralSubstitutionModelParser.java
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

import dr.evomodel.substmodel.CompositeSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

import java.util.List;
import java.util.logging.Logger;

/**
 * Parses a CompositeSubstitutionModel.
 */
public class CompositeSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String COMPOSITE_SUBSTITUTION_MODEL = "compositeSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String NORMALIZED = "normalized";

    public String getParserName() {
        return COMPOSITE_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<SubstitutionModel> substitutionModels = xo.getAllChildren(SubstitutionModel.class);

        if (substitutionModels.size() != 2 ||
                !substitutionModels.get(0).getDataType().equals(substitutionModels.get(1).getDataType())) {
            throw new XMLParseException("CompositeSubstitutionModel should take 2 substitution models of the same DataType");
        }
        int stateCount = substitutionModels.get(0).getDataType().getStateCount();
        Logger.getLogger("dr.evomodel").info("  Composite Substitution Model (" +
                substitutionModels.get(0).getId() + " + " +
                substitutionModels.get(1).getId() + ")");

        CompositeSubstitutionModel model = new CompositeSubstitutionModel(getParserName(),
                substitutionModels.get(0).getDataType(),
                substitutionModels.get(0),
                substitutionModels.get(1)
        );

//        if (!xo.getAttribute(NORMALIZED, true)) {
//            model.setNormalization(false);
//            Logger.getLogger("dr.app.beagle.evomodel").info("\tNormalization: false");
//        }

        return model;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A substitution model that is a composite of two other substitution models.";
    }

    public Class getReturnType() {
        return CompositeSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SubstitutionModel.class),
            new ElementRule(SubstitutionModel.class),
            AttributeRule.newBooleanRule(NORMALIZED, true),
    };
}
