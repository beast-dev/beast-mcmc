/*
 * RandomWalkGeneratorParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RandomWalkGenerator;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by mkarcher on 4/5/17.
 */
public class RandomWalkGeneratorParser extends AbstractXMLObjectParser {
    public static final String RANDOM_WALK_GENERATOR = "randomWalkGenerator";
    public static final String DATA = "data";
    public static final String FIRST_ELEM_PREC = "firstElementPrecision";
    public static final String PREC = "precision";
    public static final String DIM = "dimension";

    @Override
    public String getParserName() { return RANDOM_WALK_GENERATOR; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(FIRST_ELEM_PREC);
        Parameter firstElementPrecision = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PREC);
        Parameter prec = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(DATA); // May need to adapt to multiple trees, a la CoalescentLikelihoodParser
        Parameter data = (Parameter) cxo.getChild(Parameter.class);

        return new RandomWalkGenerator(data, firstElementPrecision, prec);
    }

    @Override
    public String getParserDescription() {
        return "This element generates a regular Gaussian random walk.";
    }

    @Override
    public Class getReturnType() {
        return RandomWalkGenerator.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The data to evaluate the density"),


            new ElementRule(FIRST_ELEM_PREC, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The precision for the first element of the regular random walk"),

            new ElementRule(PREC, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The precision for the relationship between adjacent elements in the random walk"),

//            AttributeRule.newIntegerRule(DIM)
    };
}
