/*
 * MicrosatelliteSamplerTreeLikelihoodParser.java
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

package dr.oldevomodelxml.treelikelihood;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.oldevomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.oldevomodel.treelikelihood.MicrosatelliteSamplerTreeLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MicrosatelliteSamplerTreeLikelihood
 */
@Deprecated // Switching to BEAGLE
public class MicrosatelliteSamplerTreeLikelihoodParser extends AbstractXMLObjectParser {
    public static final String TREE_LIKELIHOOD = "microsatelliteSamplerTreeLikelihood";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String BRANCH_RATE_MODEL = "branchRateModel";
    public String getParserName(){
        return TREE_LIKELIHOOD;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MicrosatelliteSamplerTreeModel mstModel = (MicrosatelliteSamplerTreeModel)xo.getChild(MicrosatelliteSamplerTreeModel.class);
        MicrosatelliteModel microsatelliteModel = (MicrosatelliteModel)xo.getChild(MicrosatelliteModel.class);

        BranchRateModel branchRateModel;
        Object cxo = xo.getChild(BranchRateModel.class);

        if(xo.getChild(BranchRateModel.class) !=null){

            branchRateModel = (BranchRateModel)cxo;
            System.out.println("BranchRateModel provided to MicrosatelliteSamplerTreeLikelihood");

        }else if(xo.hasChildNamed(MUTATION_RATE)){

            Parameter muRate = (Parameter)xo.getElementFirstChild(MUTATION_RATE);
            branchRateModel = new StrictClockBranchRates(muRate);
            System.out.println("mutation rate provided to MicrosatelliteSamplerTreeLikelihood");

        }else{
            branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));
        }

        return new MicrosatelliteSamplerTreeLikelihood(mstModel,microsatelliteModel, branchRateModel);
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MicrosatelliteSamplerTreeModel.class),
            new ElementRule(MicrosatelliteModel.class),
            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(BranchRateModel.class,true)
    };

    public String getParserDescription(){
        return "this parser returns an object of the TreeMicrosatelliteSamplerLikelihood class";
    }

    public Class getReturnType(){
        return MicrosatelliteSamplerTreeLikelihood.class;
    }


}
