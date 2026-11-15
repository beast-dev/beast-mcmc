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

import java.util.logging.Logger;
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
            MATDataLikelihoodDelegate delegate = getDelegate(xo);
            boolean local;
            if (xo.hasAttribute(BLOCK)) {
                String block = xo.getStringAttribute(BLOCK);
                if (!block.equalsIgnoreCase(PATH_BLOCK) && !block.equalsIgnoreCase(LOCAL_BLOCK)) {
                    throw new XMLParseException("The " + BLOCK + " attribute of " + xo.getId() + " must be '" +
                            PATH_BLOCK + "' or '" + LOCAL_BLOCK + "'");
                }
                local = block.equalsIgnoreCase(LOCAL_BLOCK);
            } else {
                // the local block where the model allows it (it is faster and mixes as well), else the path block
                local = MutationHistoryResampler.supportsLocalBlock(delegate);
            }
            if (delegate.getResampler() != null) {
                throw new XMLParseException("The MAT likelihood already has a mutation history resampler: " + xo.getId() +
                        " would be a second one (one resampler, attached to the likelihood, serves every subtree move)");
            }
            MutationHistoryResampler resampler = new MutationHistoryResampler(delegate, local);
            Logger.getLogger("dr.evomodel").info("  Mutation history resampler " + xo.getId() + ": the " +
                    (local ? LOCAL_BLOCK : PATH_BLOCK) + " block");
            return resampler;
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
                AttributeRule.newStringRule(BLOCK, true, "Which states a move resamples: 'local' (the regrafted " +
                        "parent, and the nodes whose data at a site come only from the moved subtree; the default where " +
                        "the model allows it) or 'path' (the regrafted parent and the nodes on the path between its old " +
                        "and new positions; the default for branch history models with explicit hits, per-site rates " +
                        "or unequal pattern weights, which the local block does not support). The local block needs a " +
                        "mutationHistoryGibbsOperator with this resampler, to mix the history away from the moves"),
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
