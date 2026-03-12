/*
 * CoalescentRescalerParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentRescaler;
import dr.evomodel.coalescent.demographicmodel.ConstantPopulationModel;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.demographicmodel.RescaleAwareDemographic;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

public class CoalescentRescalerParser extends AbstractXMLObjectParser {

    public static final String COALESCENT_RESCALER = "coalescentRescaler";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        if (tree == null) {
            tree = new DefaultTreeModel((Tree) xo.getChild(Tree.class));
        }
        ConstantPopulationModel constantPopulation = (ConstantPopulationModel) xo.getChild(ConstantPopulationModel.class);
        RescaleAwareDemographic rescaleAwareDemographic = (RescaleAwareDemographic) xo.getChild(RescaleAwareDemographic.class);
        CoalescentRescaler coalescentRescaler = new CoalescentRescaler(tree, constantPopulation, rescaleAwareDemographic);
        return coalescentRescaler.rescaleTree();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(new ElementRule(TreeModel.class),
                    new ElementRule(Tree.class)),
            new ElementRule(ConstantPopulationModel.class, 0, Integer.MAX_VALUE),
            new ElementRule(RescaleAwareDemographic.class, 0, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return "";
    }

    @Override
    public Class getReturnType() {
        return TreeModel.class;
    }

    @Override
    public String getParserName() {
        return COALESCENT_RESCALER;
    }
}
