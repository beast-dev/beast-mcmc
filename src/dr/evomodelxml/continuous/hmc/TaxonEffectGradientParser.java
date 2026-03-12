/*
 * TaxonEffectGradientParser.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.hmc.IntegratedLoadingsGradient;
import dr.evomodel.continuous.hmc.TaxonEffectGradient;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.inference.model.Parameter;
import dr.util.TaskPool;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.continuous.hmc.TaxonEffectGradient.Partition;

public class TaxonEffectGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "taxonEffectGradient";
    private static final String BLOCK = "partition";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<XMLObject> cxos = xo.getAllChildren(BLOCK);

        List<Partition> partitions = new ArrayList<>();
        for (XMLObject cxo : cxos) {

            TreeDataLikelihood likelihood = (TreeDataLikelihood) cxo.getChild(TreeDataLikelihood.class);
            if (!(likelihood.getDataLikelihoodDelegate() instanceof ContinuousDataLikelihoodDelegate)) {
                throw new XMLParseException("Only accepts continuous diffusion models");
            }

            partitions.add(new Partition(
                    likelihood, (ContinuousDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate(),
                    (TaxonEffectTraitDataModel) cxo.getChild(TaxonEffectTraitDataModel.class)));
        }

        Parameter effects = partitions.get(0).model.getEffects();
        for (Partition p : partitions) {
            if (p.model.getEffects() != effects) {
                throw new XMLParseException("Effects parameters must match");
            }
        }

        TaskPool taskPool = null;
        TaxonEffectGradient.ThreadUseProvider threadUse = TaxonEffectGradient.ThreadUseProvider.SERIAL;

        return new TaxonEffectGradient(partitions, taskPool, threadUse);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Generates a gradient provider for the loadings matrix when factors are integrated out";
    }

    @Override
    public Class getReturnType() {
        return IntegratedLoadingsGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    protected final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BLOCK, new XMLSyntaxRule[]{
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(TaxonEffectTraitDataModel.class),
            }, 1, Integer.MAX_VALUE),
    };
}
