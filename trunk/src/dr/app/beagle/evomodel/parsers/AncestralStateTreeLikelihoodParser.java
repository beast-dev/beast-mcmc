/*
 * AncestralStateTreeLikelihoodParser.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Map;
import java.util.Set;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */

public class AncestralStateTreeLikelihoodParser extends BeagleTreeLikelihoodParser {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
    public static final String RECONSTRUCTION_TAG = AncestralStateTreeLikelihood.STATES_KEY;
    public static final String RECONSTRUCTION_TAG_NAME = "stateTagName";
    public static final String MAP_RECONSTRUCTION = "useMAP";
    public static final String MARGINAL_LIKELIHOOD = "useMarginalLikelihood";

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

		
//		System.err.println("XML object: " + xo.toString());
	
		DataType dataType = branchModel.getRootSubstitutionModel().getDataType();

        // default tag is RECONSTRUCTION_TAG
        String tag = xo.getAttribute(RECONSTRUCTION_TAG_NAME, RECONSTRUCTION_TAG);

        boolean useMAP = xo.getAttribute(MAP_RECONSTRUCTION, false);
        boolean useMarginalLogLikelihood = xo.getAttribute(MARGINAL_LIKELIHOOD, true);

        return new AncestralStateBeagleTreeLikelihood(  // Current just returns a OldBeagleTreeLikelihood
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                partialsRestrictions,
                dataType,
                tag,
                useMAP,
                useMarginalLogLikelihood
        );
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newBooleanRule(OldTreeLikelihoodParser.USE_AMBIGUITIES, true),
                AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(GammaSiteRateModel.class),
                new ElementRule(BranchModel.class, true),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(TipStatesModel.class, true),
                new ElementRule(SubstitutionModel.class, true),
                AttributeRule.newStringRule(OldTreeLikelihoodParser.SCALING_SCHEME,true),
                new ElementRule(PARTIALS_RESTRICTION, new XMLSyntaxRule[] {
                        new ElementRule(TaxonList.class),
                        new ElementRule(Parameter.class),
                }, true),
        };
    }
}
