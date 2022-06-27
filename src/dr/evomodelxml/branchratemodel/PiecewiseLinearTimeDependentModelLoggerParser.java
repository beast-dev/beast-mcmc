/*
 * PiecewiseLinearTimeDependentModelParser.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.PiecewiseLinearTimeDependentModel;
import dr.evomodel.branchratemodel.PiecewiseLinearTimeDependentModelLogger;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class PiecewiseLinearTimeDependentModelLoggerParser extends AbstractXMLObjectParser {

    public static final String STATISTIC_NAME = "piecewiseLinearTimeEffectLogger";

    public String getParserName() {
        return STATISTIC_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        PiecewiseLinearTimeDependentModel model = (PiecewiseLinearTimeDependentModel)
                xo.getChild(PiecewiseLinearTimeDependentModel.class);

        return new PiecewiseLinearTimeDependentModelLogger(model);
    }
    
    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A parser";
    }

    public Class getReturnType() {
        return PiecewiseLinearTimeDependentModelLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PiecewiseLinearTimeDependentModel.class),
    };
}
