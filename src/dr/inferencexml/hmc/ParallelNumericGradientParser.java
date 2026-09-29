/*
 * GradientWrapperParser.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ParallelNumericalGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class ParallelNumericGradientParser extends AbstractXMLObjectParser {

    public static final String NAME = "parallelNumericGradient";
    public static final String MAIN = "main";
    public static final String WORKER = "worker";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Likelihood> likelihoods = new ArrayList<>();
        List<Parameter> parameters = new ArrayList<>();

        XMLObject main = xo.getChild(MAIN);
        likelihoods.add((Likelihood)
                main.getChild(Likelihood.class));
        parameters.add((Parameter) main.getChild(Parameter.class));

        for (XMLObject worker : xo.getAllChildren(WORKER)) {
            likelihoods.add((Likelihood)
                    worker.getChild(Likelihood.class));
            parameters.add((Parameter) worker.getChild(Parameter.class));
        }


        return new ParallelNumericalGradient(likelihoods, parameters);
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(MAIN, new XMLSyntaxRule[]{
                    new ElementRule(Likelihood.class),
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(WORKER, new XMLSyntaxRule[]{
                    new ElementRule(Likelihood.class),
                    new ElementRule(Parameter.class),
            }, 0, Integer.MAX_VALUE),
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GradientWrtParameterProvider.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}