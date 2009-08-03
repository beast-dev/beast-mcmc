/*
 * CompoundBranchRateModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.*;
import dr.xml.*;

import java.util.logging.Logger;
import java.util.Collection;
import java.util.ArrayList;
import java.util.List;

/**
 * Takes a collection of BranchRateModels and returns the product of the respective rates. In order for this to
 * work, one model should drive the actual rate of evolution and the others should be set up to provide
 * relative rates.
 * @author Andrew Rambaut
 * @version $Id:=$
 */
public class CompoundBranchRateModel extends AbstractModel implements BranchRateModel {

    public static final String COMPOUND_BRANCH_RATE_MODEL = "compoundBranchRateModel";

    private final List<BranchRateModel> branchRateModels = new ArrayList<BranchRateModel>();

    public CompoundBranchRateModel(Collection<BranchRateModel> branchRateModels) {
        super(COMPOUND_BRANCH_RATE_MODEL);
        for (BranchRateModel branchRateModel : branchRateModels) {
            addModel(branchRateModel);
            this.branchRateModels.add(branchRateModel);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(Tree tree, NodeRef node) {
        double rate = 1.0;
        for (BranchRateModel branchRateModel : branchRateModels) {
            rate *= branchRateModel.getBranchRate(tree, node);
        }
        return rate;
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPOUND_BRANCH_RATE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<BranchRateModel> branchRateModels = new ArrayList<BranchRateModel>();
            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof BranchRateModel) {
                    branchRateModels.add((BranchRateModel) xo.getChild(i));
                } else {

                    Object rogueElement = xo.getChild(i);

                    throw new XMLParseException("An element (" + rogueElement + ") which is not a branchRateModel has been added to a " + COMPOUND_BRANCH_RATE_MODEL + " element");
                }

            }

            Logger.getLogger("dr.evomodel").info("Creating a compound branch rate model of " + branchRateModels.size() + " sub-models");

            return new CompoundBranchRateModel(branchRateModels);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element provides a strict clock model. " +
                            "All branches have the same rate of molecular evolution.";
        }

        public Class getReturnType() {
            return CompoundBranchRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(BranchRateModel.class, "The component branch rate models", 1, Integer.MAX_VALUE),
        };
    };

}