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
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.List;
import java.util.logging.Logger;

/**
 * Parses a CompositeSubstitutionModel.
 */
public class CompositeSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String COMPOSITE_SUBSTITUTION_MODEL = "compositeSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String MODEL = "model";
    public static final String NORMALIZED = "normalized";

    public String getParserName() {
        return COMPOSITE_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<XMLObject> cxos = xo.getAllChildren(MODEL);
        if (cxos.size() != 2) {
            throw new XMLParseException("CompositeSubstitutionModel should take 2 GeneralSubstitutionModel of the same DataType");
        }

        GeneralSubstitutionModel substitutionModel1 = (GeneralSubstitutionModel)cxos.get(0).getChild(GeneralSubstitutionModel.class);
        Parameter weightParameter1 = (Parameter)cxos.get(0).getChild(Parameter.class);

        GeneralSubstitutionModel substitutionModel2 = (GeneralSubstitutionModel)cxos.get(1).getChild(GeneralSubstitutionModel.class);
        Parameter weightParameter2 = (Parameter)cxos.get(1).getChild(Parameter.class);
        if (!substitutionModel1.getDataType().equals(substitutionModel2.getDataType())) {
            throw new XMLParseException("GeneralSubstitutionModels are not of the same DataType");
        }

        int stateCount = substitutionModel1.getDataType().getStateCount();
        Logger.getLogger("dr.evomodel").info("  Composite Substitution Model (" +
                substitutionModel1.getId() + " + " +
                substitutionModel2.getId() + ")");

        CompositeSubstitutionModel model = new CompositeSubstitutionModel(getParserName(),
                substitutionModel1.getDataType(),
                substitutionModel1,
                weightParameter1,
                substitutionModel2,
                weightParameter2
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
            AttributeRule.newBooleanRule(NORMALIZED, true),
            new ElementRule(MODEL, new XMLSyntaxRule[] {
                    new ElementRule(GeneralSubstitutionModel.class),
                    new ElementRule(Parameter.class)
            }, 2, 2)
    };
}
