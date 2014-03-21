/*
 * LatentStateBranchRateModelParser.java
 *
 * Copyright (C) 2002-2014 Alexei Drummond, Andrew Rambaut & Marc Suchard
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.CompoundBranchRateModel;
import dr.evomodel.branchratemodel.LatentStateBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class LatentStateBranchRateModelParser extends AbstractXMLObjectParser {
    public static final String LATENT_TRANSITION_RATE = "latentTransitionRate";
    public static final String LATENT_STATE_PROPORTIONS = "latentStateProportions";


    public String getParserName() {
        return LatentStateBranchRateModel.LATENT_STATE_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        BranchRateModel nonLatentRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter latentTransitionRateParameter = (Parameter) xo.getElementFirstChild(LATENT_TRANSITION_RATE);

        Parameter latentStateProportionParameter = null;
        if (xo.hasChildNamed(LATENT_STATE_PROPORTIONS)) {
            latentStateProportionParameter = (Parameter) xo.getElementFirstChild(LATENT_STATE_PROPORTIONS);
        }

        Logger.getLogger("dr.evomodel").info("Creating a latent state branch rate model");

        return new LatentStateBranchRateModel(tree, nonLatentRateModel, latentTransitionRateParameter, latentStateProportionParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides a model with a latent state where no evolution occurs but condition on being non-latent at the nodes.";
    }

    public Class getReturnType() {
        return LatentStateBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            new ElementRule(BranchRateModel.class, "A branch rate model to provide the rates for the non-latent state"),
            new ElementRule(TreeModel.class, "The tree on which this will operate"),
            new ElementRule(LATENT_TRANSITION_RATE, Parameter.class, "A parameter which gives the instantaneous rate of switching to and from the latent state", false),
            new ElementRule(LATENT_STATE_PROPORTIONS, Parameter.class, "The proportion of each branch which is spend in a latent state", true)

    };

}
