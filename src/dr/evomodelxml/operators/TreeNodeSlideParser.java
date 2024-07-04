/*
 * TreeNodeSlideParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.TreeNodeSlide;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TreeNodeSlideParser extends AbstractXMLObjectParser {
    public static final String TREE_NODE_REHEIGHT = "nodeReHeight";

    public String getParserName() {
        return TREE_NODE_REHEIGHT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SpeciesBindings species = (SpeciesBindings) xo.getChild(SpeciesBindings.class);
        SpeciesTreeModel tree = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
//            final double range = xo.getAttribute("range", 1.0);
//            if( range <= 0 || range > 1.0 ) {
//                throw new XMLParseException("range out of range");
//            }
        //final boolean oo = xo.getAttribute("outgroup", false);
        return new TreeNodeSlide(tree, species /*, range*//*, oo*/, weight);
    }

    public String getParserDescription() {
        return "Specialized Species tree operator, transform tree without breaking embedding of gene trees.";
    }

    public Class getReturnType() {
        return TreeNodeSlide.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
               // AttributeRule.newDoubleRule("range", true),
              //  AttributeRule.newBooleanRule("outgroup", true),

                new ElementRule(SpeciesBindings.class),
                new ElementRule(SpeciesTreeModel.class)
        };
    }

}
