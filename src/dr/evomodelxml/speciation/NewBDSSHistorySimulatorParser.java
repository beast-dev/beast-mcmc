/*
 * NewBDSSHistorySimulatorParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.NewBDSSHistorySimulator;
import dr.evomodel.speciation.NewBirthDeathSerialSamplingModel;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

public class NewBDSSHistorySimulatorParser extends AbstractXMLObjectParser {

    public static final String NAME = "newBDSSHistorySimulator";
    public static final String CONDITION_SURVIVAL = "conditionOnSurvival";
    public static final String ORIGIN = "startAtOrigin";
    public static final String RECORD = "recordBeforeSampling";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        NewBirthDeathSerialSamplingModel bdss = (NewBirthDeathSerialSamplingModel) xo.getChild(NewBirthDeathSerialSamplingModel.class);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        boolean survival = xo.getAttribute(CONDITION_SURVIVAL, false);
        boolean origin = xo.getAttribute(ORIGIN, true);
        boolean recordFirst = xo.getAttribute(RECORD, true);

        return new NewBDSSHistorySimulator(bdss, tree, origin, survival, recordFirst);
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(NewBirthDeathSerialSamplingModel.class, "The birth-death model."),
            new ElementRule(Tree.class,true),
            AttributeRule.newBooleanRule("conditionOnSurvival",true,"Should only simulations which survive until the time of the most recent sample (the present) be considered? Default: false."),
            AttributeRule.newBooleanRule("startAtOrigin",true,"Should the simulation start at the origin (with one lineage, true) or at the root (with two lineages, false)? Default: true."),
            AttributeRule.newBooleanRule("recordBeforeSampling",true,"For models which have rho > 0 and r > 0, should n(t) be recorded before the intensive-sampling event hits (as opposed to after)? Default: true."),
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Simulates model-predicted counts of infected individuals at grid times. Does not simulate a complete history of every birth, death, and sample.";
    }

    @Override
    public Class<NewBDSSHistorySimulator> getReturnType() {
        return NewBDSSHistorySimulator.class;
    }

}// END: class
