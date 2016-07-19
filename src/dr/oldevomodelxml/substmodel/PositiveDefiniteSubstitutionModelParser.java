/*
 * PositiveDefiniteSubstitutionModelParser.java
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

package dr.oldevomodelxml.substmodel;

import dr.oldevomodel.substmodel.PositiveDefiniteSubstitutionModel;
import dr.inference.model.MatrixParameter;
import dr.xml.*;

/**
 *
 */
public class PositiveDefiniteSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String SVS_GENERAL_SUBSTITUTION_MODEL = "positiveDefiniteSubstitutionModel";

    public String getParserName() {
        return SVS_GENERAL_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameter ratesParameter = (MatrixParameter) xo.getChild(MatrixParameter.class);

        return new PositiveDefiniteSubstitutionModel(ratesParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of sequence substitution for any data type with stochastic variable selection.";
    }

    public Class getReturnType() {
        return PositiveDefiniteSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(MatrixParameter.class)
    };

}
