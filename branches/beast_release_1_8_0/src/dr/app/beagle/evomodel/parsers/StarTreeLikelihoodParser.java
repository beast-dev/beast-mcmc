/*
 * StarTreeLikelihoodParser.java
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

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.treelikelihood.OldBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.evomodel.treelikelihood.StarTreeLikelihood;
import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.Map;
import java.util.Set;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */

public class StarTreeLikelihoodParser extends OldTreeLikelihoodParser {

    public static final String STAR_TREE = "starTreeLikelihood";
//    public static final String RECONSTRUCTION_TAG = AncestralStateTreeLikelihood.STATES_KEY;
//    public static final String RECONSTRUCTION_TAG_NAME = "stateTagName";
//    public static final String MAP_RECONSTRUCTION = "useMAP";
//    public static final String MARGINAL_LIKELIHOOD = "useMarginalLikelihood";

    public String getParserName() {
        return STAR_TREE;
    }

    protected OldBeagleTreeLikelihood createTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                                     BranchSubstitutionModel branchSubstitutionModel, GammaSiteRateModel siteRateModel,
                                                     BranchRateModel branchRateModel,
                                                     boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                     Map<Set<String>, Parameter> partialsRestrictions,
                                                     XMLObject xo) throws XMLParseException {
           return new StarTreeLikelihood(
                    patternList,
                    treeModel,
                   branchSubstitutionModel,
                    siteRateModel,
                    branchRateModel,
                    useAmbiguities,
                    scalingScheme,
                    partialsRestrictions
            );
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return OldTreeLikelihoodParser.rules;
    }
}
