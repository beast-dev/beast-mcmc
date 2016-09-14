/*
 * ContinuousDataLikelihoodParser.java
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

package dr.evomodelxml.treedatalikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.tree.MultivariateTraitTree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import javax.swing.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class ContinuousDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String CONJUGATE_ROOT_PRIOR = AbstractMultivariateTraitLikelihood.CONJUGATE_ROOT_PRIOR;
    public static final String USE_TREE_LENGTH = AbstractMultivariateTraitLikelihood.USE_TREE_LENGTH;
    public static final String SCALE_BY_TIME = AbstractMultivariateTraitLikelihood.SCALE_BY_TIME;
    public static final String RECIPROCAL_RATES = AbstractMultivariateTraitLikelihood.RECIPROCAL_RATES;
    public static final String PRIOR_SAMPLE_SIZE = AbstractMultivariateTraitLikelihood.PRIOR_SAMPLE_SIZE;

    public static final String CONTINUOUS_DATA_LIKELIHOOD = "traitDataLikelihood";

    public String getParserName() {
        return CONTINUOUS_DATA_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
        String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, true);
        CompoundParameter traitParameter = returnValue.traitParameter;
        List<Integer> missingIndices = returnValue.missingIndices;
        traitName = returnValue.traitName;

        ContinuousTraitDataModel dataModel = new ContinuousTraitDataModel(traitName,
                traitParameter,
                missingIndices);

        ContinuousDataLikelihoodDelegate delegate = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionModel, dataModel, rateModel);

        return new TreeDataLikelihood(delegate, treeModel, rateModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of trait data on a tree given a diffusion model.";
    }

    public Class getReturnType() {
        return TreeDataLikelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(MultivariateDiffusionModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(CONJUGATE_ROOT_PRIOR, ConjugateRootTraitPrior.rules),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
