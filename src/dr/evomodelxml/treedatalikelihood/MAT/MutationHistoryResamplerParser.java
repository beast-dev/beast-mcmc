/*
 * MutationHistoryResamplerParser.java
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

import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MAT.MATDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.xml.*;
import dr.evomodel.treedatalikelihood.MAT.MutationHistoryResampler;

/**
 * Parses a mutation history resampler, which can be given to subtreeLeap, subtreeJump and
 * fixedHeightSubtreePruneRegraft operators to resample the MAT along with their topology moves.
 */
public class MutationHistoryResamplerParser extends AbstractXMLObjectParser {

    public static final String MUTATION_HISTORY_RESAMPLER = "mutationHistoryResampler";

    @Override
    public String getParserName() {
        return MUTATION_HISTORY_RESAMPLER;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new MutationHistoryResampler(getDelegate(xo));
    }

    static MATDataLikelihoodDelegate getDelegate(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        DataLikelihoodDelegate delegate = likelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof MATDataLikelihoodDelegate)) {
            throw new XMLParseException("Element " + xo.getId() + " requires a " +
                    MATDataLikelihoodParser.MAT_DATA_LIKELIHOOD + " but " + likelihood.getId() + " is not one");
        }
        return (MATDataLikelihoodDelegate) delegate;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(TreeDataLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Resamples the ancestral states of a mutation annotated tree jointly with subtree prune-and-regraft moves.";
    }

    @Override
    @SuppressWarnings("rawtypes")
    public Class getReturnType() {
        return MutationHistoryResampler.class;
    }
}
