/*
 * HierarchicalSubstitutionModelParser.java
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Hierarchical Rate Substitution Model
 *
 * @author Jeff Thorne
 * @author Tae-Kun Seo
 * @author Sangchul Choi
 * @author Xiang Ji
 */
public class HierarchicalSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String HIERARCHICAL_SUBSTITUTION_MODEL = "hierarchicalSubstitutionModel";
    public static final String CHARACTER_TREE = "characterTree";
    public static final String RATES_PARAMETER = "ratesParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        TreeModel tree = (TreeModel) xo.getElementFirstChild(CHARACTER_TREE);

        HierarchicalSubstitutionModel.HierarchicalRateProvider hierarchicalRateProvider = xo.hasChildNamed(RATES_PARAMETER) ?
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default((Parameter) xo.getElementFirstChild(RATES_PARAMETER), tree)
                : new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(new Parameter.Default(tree.getNodeCount(), 1.0), tree);

        return new HierarchicalSubstitutionModel(HIERARCHICAL_SUBSTITUTION_MODEL, hierarchicalRateProvider, freqModel);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FrequencyModelParser.FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(CHARACTER_TREE,
                    new XMLSyntaxRule[]{new ElementRule(TreeModel.class)}),
            new ElementRule(RATES_PARAMETER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };

    @Override
    public String getParserDescription() {
        return "";
    }

    @Override
    public Class getReturnType() {
        return HierarchicalSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return HIERARCHICAL_SUBSTITUTION_MODEL;
    }
}
