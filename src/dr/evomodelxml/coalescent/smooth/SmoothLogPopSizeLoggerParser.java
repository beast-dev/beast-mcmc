/*
 * SmoothLogPopSizeLoggerParser.java
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

package dr.evomodelxml.coalescent.smooth;

import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.coalescent.smooth.NewSmoothSkygridLikelihood;
import dr.evomodel.coalescent.smooth.SmoothLogPopSizeLogger;
import dr.xml.*;

public class SmoothLogPopSizeLoggerParser extends AbstractXMLObjectParser {

    private static final String SMOOTH_LOG_POP_SIZE_PARAMETER = "smoothLogPopSizeParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Object likelihood = xo.getChild(AbstractCoalescentLikelihood.class);
        // <smoothSkygridLikelihood> builds a NewSmoothSkygridLikelihood only with implementation="new"
        if (!(likelihood instanceof NewSmoothSkygridLikelihood)) {
            throw new XMLParseException("This element needs a <smoothSkygridLikelihood> with implementation=\"new\"");
        }
        return new SmoothLogPopSizeLogger((NewSmoothSkygridLikelihood) likelihood);
    }

    @Override
    public String getParserName() {
        return SMOOTH_LOG_POP_SIZE_PARAMETER;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(AbstractCoalescentLikelihood.class)
    };

    @Override
    public String getParserDescription() {
        return "Logs the log of smoothed pop sizes of smooth Skygrid model";
    }

    @Override
    public Class getReturnType() {
        return SmoothLogPopSizeLogger.class;
    }
}
