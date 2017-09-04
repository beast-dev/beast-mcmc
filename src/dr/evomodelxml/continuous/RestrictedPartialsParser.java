/*
 * RestrictedPartialsParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.RestrictedPartials;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.MonophylyStatisticParser;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class RestrictedPartialsParser extends AbstractXMLObjectParser {

    public static final String RESTRICTED_PARTIALS = "restrictedPartials";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getId();

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        TaxonList taxa = MonophylyStatisticParser.parseTaxonListOrTaxa(xo.getChild(MonophylyStatisticParser.MRCA));

        Parameter meanParameter = (Parameter) xo.getElementFirstChild(MultivariateDistributionLikelihood.MVN_MEAN);
        Parameter priorSampleSize =
                (Parameter) xo.getElementFirstChild(AbstractMultivariateTraitLikelihood.PRIOR_SAMPLE_SIZE);

        RestrictedPartials rp = null;

        try {
            rp =  new RestrictedPartials(name, tree, taxa, meanParameter, priorSampleSize);
        } catch (TreeUtils.MissingTaxonException e) {
            throw new XMLParseException("Unable to find taxa for " + xo.getId());
        }

        return rp;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(MonophylyStatisticParser.MRCA, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }),
            new ElementRule(MultivariateDistributionLikelihood.MVN_MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(AbstractMultivariateTraitLikelihood.PRIOR_SAMPLE_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    public String getParserDescription() {
        return "A restricted partials assignment for a set of taxa";
    }

    public Class getReturnType() {
        return RestrictedPartials.class;
    }

    public String getParserName() {
        return RESTRICTED_PARTIALS;
    }

}
