/*
 * LKJCorrelationWithStructuralZerosDistributionParser.java
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

import dr.inference.model.ParameterParser;
import dr.math.distributions.LKJCorrelationWithStructuralZerosDistribution;
import dr.xml.*;

import java.util.ArrayList;

public class LKJCorrelationWithStructuralZerosDistributionParser extends AbstractXMLObjectParser {

    private static final String BLOCKS = "blocks";
    private static final String BLOCK = "block"; //name doesn't actually matter

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject bxo = xo.getChild(BLOCKS);
        ArrayList<int[]> blocks = new ArrayList<>();

        for (int i = 0; i < bxo.getChildCount(); i++) {
            XMLObject bcxo = (XMLObject) bxo.getChild(i);
            blocks.add(bcxo.getIntegerArrayChild(0));
        }

        for (int[] block : blocks) {
            for (int i = 0; i < block.length; i++) {
                block[i] = block[i] - 1;
            }
        }

        int dim = xo.getIntegerAttribute(ParameterParser.DIMENSION);
        double shape = xo.getDoubleAttribute(PriorParsers.SHAPE);

        return new LKJCorrelationWithStructuralZerosDistribution(dim, shape, blocks);

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(PriorParsers.SHAPE),
                AttributeRule.newIntegerRule(ParameterParser.DIMENSION),
                new ElementRule(BLOCKS, new XMLSyntaxRule[]{
                        new ElementRule(BLOCK, new XMLSyntaxRule[0], "", 0, Integer.MAX_VALUE)
                })

        };
    }

    @Override
    public String getParserDescription() {
        return "LKJ correlation distribution with some diagonal blocks fixed at the identity matrix";
    }

    @Override
    public Class getReturnType() {
        return LKJCorrelationWithStructuralZerosDistribution.class;
    }

    @Override
    public String getParserName() {
        return LKJCorrelationWithStructuralZerosDistribution.LKJ_WITH_ZEROS;
    }
}
