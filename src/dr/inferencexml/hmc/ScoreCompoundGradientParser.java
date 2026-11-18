/*
 * ScoreCompoundGradientParser.java
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

package dr.inferencexml.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ScoreCompoundGradient;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Parses {@code <scoreCompoundGradient>}: a compound gradient over its child gradients that also
 * concatenates their score contributions (see {@link ScoreCompoundGradient}).  Use it in place of
 * {@code <compoundGradient>} where an empirical Fisher information is wanted, for example as the
 * {@code <fisherGradient>} of {@code <maximizeWrtParameter>}.
 */
public class ScoreCompoundGradientParser extends AbstractXMLObjectParser {

    public static final String SCORE_COMPOUND_GRADIENT = "scoreCompoundGradient";

    @Override
    public String getParserName() {
        return SCORE_COMPOUND_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final List<GradientWrtParameterProvider> gradients = new ArrayList<GradientWrtParameterProvider>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            gradients.add((GradientWrtParameterProvider) xo.getChild(i));
        }
        return new ScoreCompoundGradient(gradients);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return "A compound gradient that also exposes concatenated score contributions.";
    }

    @Override
    public Class getReturnType() {
        return ScoreCompoundGradient.class;
    }
}
