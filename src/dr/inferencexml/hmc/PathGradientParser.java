/*
 * JointGradientParser.java
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

package dr.inferencexml.hmc;

import dr.inference.hmc.PathGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.PathLikelihood;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class PathGradientParser extends AbstractXMLObjectParser {

    private final static String PATH_GRADIENT = "pathGradient";
    
    @Override
    public String getParserName() {
        return PATH_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GradientWrtParameterProvider source = (GradientWrtParameterProvider)
                xo.getElementFirstChild(PathLikelihood.SOURCE);

        GradientWrtParameterProvider destination = (GradientWrtParameterProvider)
                xo.getElementFirstChild(PathLikelihood.DESTINATION);

        return new PathGradient(source, destination);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PathLikelihood.SOURCE, new XMLSyntaxRule[] {
                    new ElementRule(GradientWrtParameterProvider.class),
            }),
            new ElementRule(PathLikelihood.DESTINATION, new XMLSyntaxRule[] {
                    new ElementRule(GradientWrtParameterProvider.class),
            }),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PathGradient.class;
    }
}
