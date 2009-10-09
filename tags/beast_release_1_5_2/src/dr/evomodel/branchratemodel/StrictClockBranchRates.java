/*
 * StrictClockBranchRates.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: StrictClockBranchRates.java,v 1.3 2006/01/09 17:44:30 rambaut Exp $
 */
public class StrictClockBranchRates extends AbstractModel implements BranchRateModel {

    public static final String STRICT_CLOCK_BRANCH_RATES = "strictClockBranchRates";
    public static final String RATE = "rate";

    private final Parameter rateParameter;

    public StrictClockBranchRates(Parameter rateParameter) {

        super(STRICT_CLOCK_BRANCH_RATES);

        this.rateParameter = rateParameter;

        addVariable(rateParameter);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // nothing to do
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
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
        return rateParameter.getParameterValue(0);
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return STRICT_CLOCK_BRANCH_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);

            Logger.getLogger("dr.evomodel").info("Using strict molecular clock model.");

            return new StrictClockBranchRates(rateParameter);
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
            return StrictClockBranchRates.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(RATE, Parameter.class, "The molecular evolutionary rate parameter", false),
        };
    };

}