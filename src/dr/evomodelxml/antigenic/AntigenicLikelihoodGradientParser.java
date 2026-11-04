/*
 * AntigenicLikelihoodGradientParser.java
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

package dr.evomodelxml.antigenic;

import dr.evomodel.antigenic.AntigenicGradientWrtParameter;
import dr.evomodel.antigenic.AntigenicLikelihoodGradient;
import dr.evomodel.antigenic.NewAntigenicLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class AntigenicLikelihoodGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "antigenicGradient";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        NewAntigenicLikelihood likelihood = (NewAntigenicLikelihood) xo.getChild(NewAntigenicLikelihood.class);
        List<Parameter> parameters = xo.getAllChildren(Parameter.class);

        List<AntigenicGradientWrtParameter> wrtList = new ArrayList<>();
        for (Parameter parameter : parameters) {
            AntigenicGradientWrtParameter wrt = likelihood.wrtFactory(parameter);
            wrtList.add(wrt);
        }

        return new AntigenicLikelihoodGradient(likelihood, wrtList);
    }

    public String getParserDescription() {
        return "Provides the gradient of the likelihood of immunological assay data such as " +
                "hemagglutinin inhibition (HI) given vectors of coordinates " +
                "for viruses and sera/antisera in some multidimensional 'antigenic' space.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(NewAntigenicLikelihood.class),
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
    };

    public Class getReturnType() {
        return AntigenicLikelihoodGradient.class;
    }
}

