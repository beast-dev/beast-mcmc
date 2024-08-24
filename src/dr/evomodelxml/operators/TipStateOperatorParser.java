/*
 * TipStateOperatorParser.java
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
 */

package dr.evomodelxml.operators;

import dr.evomodel.operators.TipStateOperator;
import dr.evomodel.tipstatesmodel.TimeVaryingFrequenciesModel;
import dr.evomodel.treedatalikelihood.TipStateAccessor;
import dr.xml.*;

import java.util.List;

import static dr.inference.operators.MCMCOperator.WEIGHT;

/**
 * @author Marc A. Suchard
 */

public class TipStateOperatorParser extends AbstractXMLObjectParser {

    public static final String TIP_STATE_OPERATOR = "tipStateOperator";

    public String getParserName() { return TIP_STATE_OPERATOR; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<TipStateAccessor> accessors = xo.getAllChildren(TipStateAccessor.class);
        TimeVaryingFrequenciesModel frequencies = (TimeVaryingFrequenciesModel)
                xo.getChild(TimeVaryingFrequenciesModel.class);
        double weight = xo.getDoubleAttribute(WEIGHT);

        return new TipStateOperator(accessors, frequencies, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator to swap tip states between two random tips.";
    }

    public Class getReturnType() {
        return TipStateOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(TipStateAccessor.class, 1, Integer.MAX_VALUE),
            new ElementRule(TimeVaryingFrequenciesModel.class),
    };
}
