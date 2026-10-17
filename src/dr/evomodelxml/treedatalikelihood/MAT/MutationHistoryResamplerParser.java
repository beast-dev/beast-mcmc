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
    public static final String SITE_SELECTION = "siteSelection";
    public static final String BLOCK = "block";
    public static final String PATH_BLOCK = "path";
    public static final String LOCAL_BLOCK = "local";

    @Override
    public String getParserName() {
        return MUTATION_HISTORY_RESAMPLER;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        try {
            if (xo.hasAttribute(SITE_SELECTION) && !xo.getStringAttribute(SITE_SELECTION).equalsIgnoreCase("all")) {
                // the hybrid selection (only sites varying around the block) could not create back-mutations
                // at sites that looked uniform around every block, so it was not irreducible
                throw new XMLParseException("The " + SITE_SELECTION + " attribute of " + xo.getId() + " must be " +
                        "'all' (every site is resampled; 'hybrid' is no longer supported)");
            }
            String block = xo.getAttribute(BLOCK, PATH_BLOCK);
            if (!block.equalsIgnoreCase(PATH_BLOCK) && !block.equalsIgnoreCase(LOCAL_BLOCK)) {
                throw new XMLParseException("The " + BLOCK + " attribute of " + xo.getId() + " must be '" +
                        PATH_BLOCK + "' or '" + LOCAL_BLOCK + "'");
            }
            return new MutationHistoryResampler(getDelegate(xo), block.equalsIgnoreCase(LOCAL_BLOCK));
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage() + " in " + xo.getId());
        }
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
                AttributeRule.newStringRule(SITE_SELECTION, true, "Which sites a move resamples: only 'all' " +
                        "(kept for older files)"),
                AttributeRule.newStringRule(BLOCK, true, "Which states a move resamples: 'path' (the regrafted " +
                        "parent and the nodes on the path between its old and new positions; the default) or 'local' " +
                        "(the regrafted parent, and the nodes whose data at a site come only from the moved subtree)"),
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
