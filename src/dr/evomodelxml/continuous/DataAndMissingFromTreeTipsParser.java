/*
 * DataAndMissingFromTreeTipsParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Gabe Hassler
 * @author Marc Suchard
 */

public class DataAndMissingFromTreeTipsParser extends AbstractXMLObjectParser{

        public final static String DATA_FROM_TREE_TIPS = "dataAndMissingFromTreeTips";
        public final static String DATA = "data";
        public static final String CONTINUOUS = "continuous";


        public String getParserName() {
            return DATA_FROM_TREE_TIPS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, treeModel, true);

            return returnValue;
        }

        public static final XMLSyntaxRule[] rules = {
                new ElementRule(MutableTreeModel.class),
                AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                })
        };

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Returns the data and missing indices from the tips of a tree.";
        }

        @Override
        public Class getReturnType() {
            return TreeTraitParserUtilities.TraitsAndMissingIndices.class;
        }


}