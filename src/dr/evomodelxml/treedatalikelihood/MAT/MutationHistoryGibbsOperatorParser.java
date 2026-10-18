/*
 * MutationHistoryGibbsOperatorParser.java
 *
 * Copyright © 2026, the BEAST Development Team.
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

package dr.evomodelxml.treedatalikelihood.MAT;

import dr.evomodel.treedatalikelihood.MAT.MutationHistoryGibbsOperator;
import dr.evomodel.treedatalikelihood.MAT.MutationHistoryResampler;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Parser for the Gibbs operator on a MAT's mutational history.
 */
public class MutationHistoryGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String MUTATION_HISTORY_GIBBS_OPERATOR = "mutationHistoryGibbsOperator";
    public static final String SIZE = "size";

    @Override
    public String getParserName() {
        return MUTATION_HISTORY_GIBBS_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MutationHistoryResampler resampler = (MutationHistoryResampler) xo.getChild(MutationHistoryResampler.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int size = xo.getAttribute(SIZE, 3);
        try {
            return new MutationHistoryGibbsOperator(resampler, size, weight);
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage() + " in " + xo.getId());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newIntegerRule(SIZE, true, "The largest number of nodes resampled together: a random " +
                        "internal node and up to size - 1 of its ancestors (default 3)"),
                new ElementRule(MutationHistoryResampler.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Gibbs update of the mutational history of a mutation annotated tree, without changing the tree: " +
                "the states of a random internal node and some of its ancestors are drawn from their full conditional.";
    }

    @Override
    @SuppressWarnings("rawtypes")
    public Class getReturnType() {
        return MutationHistoryGibbsOperator.class;
    }
}
