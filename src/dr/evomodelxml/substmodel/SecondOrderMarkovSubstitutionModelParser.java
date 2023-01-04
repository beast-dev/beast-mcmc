/*
 * SecondOrderMarkovSubstitutionModelParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.SecondOrderMarkovSubstitutionModel;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class SecondOrderMarkovSubstitutionModelParser extends AbstractXMLObjectParser {

    private final String NAME = "secondOrderMarkovSubstitutionModel";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        BaseSubstitutionModel baseSubstitutionModel = (BaseSubstitutionModel) xo.getChild(BaseSubstitutionModel.class);
        SecondOrderMarkovSubstitutionModel.SecondOrderMarkovPairedDataType dataType = new SecondOrderMarkovSubstitutionModel.SecondOrderMarkovPairedDataType(baseSubstitutionModel.getDataType());
        return new SecondOrderMarkovSubstitutionModel(NAME + "(" + baseSubstitutionModel.getModelName() + ")",
                dataType, baseSubstitutionModel);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BaseSubstitutionModel.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SecondOrderMarkovSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}