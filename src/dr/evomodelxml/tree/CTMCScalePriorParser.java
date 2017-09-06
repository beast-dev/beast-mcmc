/*
 * CTMCScalePriorParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.util.TaxonList;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.CTMCScalePrior;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 *
 */
public class CTMCScalePriorParser extends AbstractXMLObjectParser {
    public static final String MODEL_NAME = "ctmcScalePrior";
    public static final String SCALEPARAMETER = "ctmcScale";
    public static final String RECIPROCAL = "reciprocal";
    public static final String TRIAL = "trial";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        TaxonList taxa = (TaxonList) xo.getChild(TaxonList.class);

        Parameter ctmcScale = (Parameter) xo.getElementFirstChild(SCALEPARAMETER);
        boolean reciprocal = xo.getAttribute(RECIPROCAL, false);
        boolean trial = xo.getAttribute(TRIAL, false);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        Logger.getLogger("dr.evolution").info("Creating CTMC Scale Reference Prior model.");
        if (taxa != null) {
            Logger.getLogger("dr.evolution").info("Acting on subtree of size " + taxa.getTaxonCount());
        }
        return new CTMCScalePrior(MODEL_NAME, ctmcScale, treeModel, taxa, reciprocal, substitutionModel, trial);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the prior for CTMC scale parameter.";
    }

    public Class getReturnType() {
        return CTMCScalePrior.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(TaxonList.class, true),
            new ElementRule(SCALEPARAMETER, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SCALEPARAMETER, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newBooleanRule(RECIPROCAL, true),
            new ElementRule(SubstitutionModel.class, true),
            AttributeRule.newBooleanRule(TRIAL, true),
    };
}
