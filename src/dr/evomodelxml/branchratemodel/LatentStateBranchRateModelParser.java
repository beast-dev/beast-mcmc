/*
 * LatentStateBranchRateModelParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.branchratemodel.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class LatentStateBranchRateModelParser extends AbstractXMLObjectParser {
    public static final String LATENT_TRANSITION_RATE = "latentTransitionRate";
    public static final String LATENT_TRANSITION_FREQUENCY = "latentTransitionFrequency";
    public static final String LATENT_STATE_PROPORTIONS = "latentStateProportions";
    public static final String NON_LATENT_RATE = "nonLatentRate";
    public static final String MAXIMUM_LATENT_PERIODS="maximumNumberOfLatentPeriods";
    public static final String SCALE_TRANSITION_RATE_BY_ROOT_HEIGHT= "scaleTransitionRateByRootHeight";


    public String getParserName() {
        return LatentStateBranchRateModel.LATENT_STATE_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        BranchRateModel nonLatentRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter latentTransitionRateParameter = (Parameter) xo.getElementFirstChild(LATENT_TRANSITION_RATE);
        Parameter latentTransitionFrequencyParameter = (Parameter) xo.getElementFirstChild(LATENT_TRANSITION_FREQUENCY);

        CountableBranchCategoryProvider branchCategoryProvider = (CountableBranchCategoryProvider)xo.getChild(CountableBranchCategoryProvider.class);

        Parameter latentStateProportionParameter = null;
        if (xo.hasChildNamed(LATENT_STATE_PROPORTIONS)) {
            latentStateProportionParameter = (Parameter) xo.getElementFirstChild(LATENT_STATE_PROPORTIONS);
        }
        Parameter nonLatentRate = null;
        if(xo.hasChildNamed(NON_LATENT_RATE)){
            nonLatentRate = (Parameter) xo.getElementFirstChild(NON_LATENT_RATE);
        }

        int maximumNumberOfLatentPeriods = xo.getAttribute(MAXIMUM_LATENT_PERIODS,5);
        boolean scaleByRootHeight = xo.getAttribute(SCALE_TRANSITION_RATE_BY_ROOT_HEIGHT,false);

        Logger.getLogger("dr").info("\nCreating a latent state branch rate model");

        return new LatentStateBranchRateModel(LatentStateBranchRateModel.LATENT_STATE_BRANCH_RATE_MODEL,
                tree, nonLatentRateModel,
                latentTransitionRateParameter, latentTransitionFrequencyParameter, /* 0/1 CTMC have two parameters */
                latentStateProportionParameter, nonLatentRate ,branchCategoryProvider,maximumNumberOfLatentPeriods,scaleByRootHeight);


//        return new SericolaLatentStateBranchRateModel(SericolaLatentStateBranchRateModel.LATENT_STATE_BRANCH_RATE_MODEL,
//                tree, nonLatentRateModel,
//                latentTransitionRateParameter, latentTransitionFrequencyParameter, /* 0/1 CTMC have two parameters */
//                latentStateProportionParameter, branchCategoryProvider);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides a model with a latent state where no evolution occurs but condition on being non-latent at the nodes.";
    }

    public Class getReturnType()  {
        return SericolaLatentStateBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            new ElementRule(BranchRateModel.class, "A branch rate model to provide the rates for the non-latent state"),
            new ElementRule(TreeModel.class, "The tree on which this will operate"),
            new ElementRule(CountableBranchCategoryProvider.class, true),
            new ElementRule(LATENT_TRANSITION_RATE, Parameter.class, "A parameter which gives the instantaneous rate of switching to and from the latent state", false),
            new ElementRule(LATENT_TRANSITION_FREQUENCY, Parameter.class, "A parameter which gives the rate bias of switching to and from the latent state", false),
            new ElementRule(LATENT_STATE_PROPORTIONS, Parameter.class, "The proportion of each branch which is spend in a latent state", true),
            new ElementRule(NON_LATENT_RATE,Parameter.class, "An optional rate that will override the branch rates from the branch rate model when proportion latent is greater than 0", true),
            new ElementRule(MAXIMUM_LATENT_PERIODS,int.class,"An optional integer specifying the number of latent periods used in calculating the likelihood of the latent proportion",true),
            new ElementRule(SCALE_TRANSITION_RATE_BY_ROOT_HEIGHT,boolean.class,"An optional boolean specifying that the transition rate parameter should be normalized by the root height.",true),

    };

}
