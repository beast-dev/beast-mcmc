/*
 * ALSTreeLikelihoodParser.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.MutationDeathModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.ALSBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.MSSD.AbstractObservationProcess;
import dr.evomodel.MSSD.AnyTipObservationProcess;
import dr.evomodel.MSSD.SingleTipObservationProcess;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */

public class ALSTreeLikelihoodParser extends BeagleTreeLikelihoodParser {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.LIKE_NAME;
    public static final String RECONSTRUCTION_TAG_NAME = AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME;
    public static final String INTEGRATE_GAIN_RATE = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.INTEGRATE_GAIN_RATE;
    public static final String OBSERVATION_TYPE = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.OBSERVATION_TYPE;
    public static final String OBSERVATION_PROCESS = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.OBSERVATION_PROCESS;
    public static final String OBSERVATION_TAXON = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.OBSERVATION_TAXON;
    public static final String IMMIGRATION_RATE = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.IMMIGRATION_RATE;
    public static final String ANY_TIP = dr.evomodelxml.MSSD.ALSTreeLikelihoodParser.ANY_TIP;

    public String getParserName() {
        return RECONSTRUCTING_TREE_LIKELIHOOD;
    }

    protected BeagleTreeLikelihood createTreeLikelihood(
            PatternList patternList, //
            TreeModel treeModel, //
            BranchModel branchModel, //
            GammaSiteRateModel siteRateModel, //
            BranchRateModel branchRateModel, //
            TipStatesModel tipStatesModel, //
            boolean useAmbiguities, //
            PartialsRescalingScheme scalingScheme, //
            Map<Set<String>, //
                    Parameter> partialsRestrictions, //
            XMLObject xo //
    ) throws XMLParseException {

        boolean integrateGainRate = xo.getBooleanAttribute(INTEGRATE_GAIN_RATE);

        useAmbiguities = true; // TODO No effect

        if (scalingScheme != PartialsRescalingScheme.NONE) {
            throw new XMLParseException("No rescaling scheme is currently support by the mutation-death model " + xo.getId());
        }

        Parameter mu = ((MutationDeathModel) siteRateModel.getSubstitutionModel()).getDeathParameter();
        Parameter lam;
        if (!integrateGainRate) {
            lam = (Parameter) xo.getElementFirstChild(IMMIGRATION_RATE);
        } else {
            lam = new Parameter.Default("gainRate", 1.0, 0.001, 1.999);
        }
        AbstractObservationProcess observationProcess = null;

        Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating a BEAGLE ALSTreeLikelihood model.");
        for (int i = 0; i < xo.getChildCount(); ++i) {
            Object cxo = xo.getChild(i);
            if (cxo instanceof XMLObject && ((XMLObject) cxo).getName().equals(OBSERVATION_PROCESS)) {
                if (((XMLObject) cxo).getStringAttribute(OBSERVATION_TYPE).equals("singleTip")) {
                    String taxonName = ((XMLObject) cxo).getStringAttribute(OBSERVATION_TAXON);
                    Taxon taxon = treeModel.getTaxon(treeModel.getTaxonIndex(taxonName));
                    observationProcess = new SingleTipObservationProcess(treeModel, patternList, siteRateModel,
                            branchRateModel, mu, lam, taxon);
                    Logger.getLogger("dr.evolution").info("All traits are assumed extant in " + taxonName);
                } else {  // "anyTip" observation process
                    observationProcess = new AnyTipObservationProcess(ANY_TIP, treeModel, patternList,
                            siteRateModel, branchRateModel, mu, lam);
                    Logger.getLogger("dr.evolution").info("Observed traits are assumed to be extant in at least one tip node.");
                }

                observationProcess.setIntegrateGainRate(integrateGainRate);
            }
        }
        Logger.getLogger("dr.evolution").info("\tIf you publish results using Acquisition-Loss-Mutation (ALS) Model likelihood, please reference Alekseyenko, Lee and Suchard (2008) Syst. Biol 57: 772-784.\n---------------------------------\n");

        return new ALSBeagleTreeLikelihood(
                observationProcess,
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                partialsRestrictions
        );
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(BeagleTreeLikelihoodParser.USE_AMBIGUITIES, true),
                AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),

                AttributeRule.newBooleanRule(INTEGRATE_GAIN_RATE),
                new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),

                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(GammaSiteRateModel.class),
                new ElementRule(BranchModel.class, true),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(TipStatesModel.class, true),
                new ElementRule(SubstitutionModel.class, true),
                AttributeRule.newStringRule(BeagleTreeLikelihoodParser.SCALING_SCHEME, true),
                new ElementRule(PARTIALS_RESTRICTION, new XMLSyntaxRule[]{
                        new ElementRule(TaxonList.class),
                        new ElementRule(Parameter.class),
                }, true),

                new ElementRule(OBSERVATION_PROCESS,
                        new XMLSyntaxRule[]{AttributeRule.newStringRule(OBSERVATION_TYPE, false),
                                AttributeRule.newStringRule(OBSERVATION_TAXON, true)})

        };
    }
}
