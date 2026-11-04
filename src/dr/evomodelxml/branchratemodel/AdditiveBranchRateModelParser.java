/*
 * AdditiveBranchRateModelParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.AdditiveBranchRateModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Karthik Gangavarapu
 * @author Marc Suchard
 */
public class AdditiveBranchRateModelParser extends AbstractXMLObjectParser {

    public static String ADDITIVE_BRANCH_RATES = "additiveBranchRates";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<AbstractBranchRateModel> branchRateModels = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            branchRateModels.add((AbstractBranchRateModel) xo.getChild(i));
        }
        return new AdditiveBranchRateModel(branchRateModels);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AbstractBranchRateModel.class, 2, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return AdditiveBranchRateModel.class;
    }

    @Override
    public String getParserName() {
        return ADDITIVE_BRANCH_RATES;
    }
}
