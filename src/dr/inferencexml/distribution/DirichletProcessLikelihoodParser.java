/*
 * DirichletProcessLikelihoodParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.DirichletProcessLikelihood;
import dr.inference.model.*;
import dr.xml.*;

/**
 *
 */
public class DirichletProcessLikelihoodParser extends AbstractXMLObjectParser {

    public static final String ETA = "eta";
    public static final String CHI = "chi";

    public String getParserName() {
        return DirichletProcessLikelihood.DIRICHLET_PROCESS_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(ETA);
        Statistic etaParam = (Statistic) cxo.getChild(Statistic.class);

        cxo = xo.getChild(CHI);
        Parameter chiParameter = (Parameter) cxo.getChild(Parameter.class);

        return new DirichletProcessLikelihood(etaParam, chiParameter);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(ETA,
                    new XMLSyntaxRule[]{new ElementRule(Statistic.class)}, "Counts of N items distributed amongst K classes"),
            new ElementRule(CHI,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, "Aggregation parameter"),
    };

    public String getParserDescription() {
        return "Calculates the likelihood of some items distributed into a number of classes under a Dirichlet drocess.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }
}
