/*
 * AnyTipObservationProcessParser.java
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

package dr.oldevomodelxml.MSSD;

import dr.evolution.alignment.PatternList;
import dr.oldevomodel.MSSD.AnyTipObservationProcess;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 *
 */
public class AnyTipObservationProcessParser extends AbstractXMLObjectParser {
    public static final String MODEL_NAME = "anyTipObservationProcess";
    final static String DEATH_RATE = "deathRate";
    final static String IMMIGRATION_RATE = "immigrationRate";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter mu = (Parameter) xo.getElementFirstChild(DEATH_RATE);
        Parameter lam = (Parameter) xo.getElementFirstChild(IMMIGRATION_RATE);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        Logger.getLogger("dr.evomodel.MSSD").info("Creating AnyTipObservationProcess model. Observed traits are assumed to be extant in at least one tip node. Initial mu = " + mu.getParameterValue(0) + " initial lam = " + lam.getParameterValue(0));

        return new AnyTipObservationProcess(MODEL_NAME, treeModel, patterns, siteModel, branchRateModel, mu, lam);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the AnyTipObservationProcess for ALSTreeLikelihood calculations";
    }

    public Class getReturnType() {
        return AnyTipObservationProcess.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(PatternList.class),
            new ElementRule(SiteModel.class),
            new ElementRule(BranchRateModel.class),
            new ElementRule(DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };

}
