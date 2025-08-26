/*
 * PolymorphismAwareSubstitutionModelParser.java
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

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PolymorphismAwareDataType;
import dr.evolution.util.Taxa;
import dr.evomodel.substmodel.PolymorphismAwareSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Nicola De Maio
 * @author Ben Redelings
 * @author Marc A. Suchard
 */
public class PolymorphismAwareSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String NAME = "polymorphismAwareSubstitutionModel";
    public static final String VIRTUAL_POP_SIZE = "virtualPopSize";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SubstitutionModel baseSubstitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        int virtualPopSize = xo.getIntegerAttribute(VIRTUAL_POP_SIZE);
        SimpleAlignment alignment = (SimpleAlignment) xo.getChild(SimpleAlignment.class);
        DataType baseDataType = alignment.getDataType();
        PolymorphismAwareDataType dataType = (PolymorphismAwareDataType) DataType.getRegisteredDataTypeByName(PolymorphismAwareDataType.getDataTypeDescription(baseDataType, virtualPopSize));

        return new PolymorphismAwareSubstitutionModel(baseSubstitutionModel, dataType);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(SubstitutionModel.class),
            new ElementRule(SimpleAlignment.class),
            AttributeRule.newIntegerRule(VIRTUAL_POP_SIZE)
    };

    @Override
    public String getParserDescription() {
        return "";
    }

    @Override
    public Class getReturnType() {
        return PolymorphismAwareSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
