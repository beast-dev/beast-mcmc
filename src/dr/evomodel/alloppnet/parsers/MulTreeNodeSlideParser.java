
/*
 * MulTreeNodeSlideParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.parsers;

import dr.evomodel.alloppnet.operators.MulTreeNodeSlide;
import dr.evomodel.alloppnet.speciation.MulSpeciesBindings;
import dr.evomodel.alloppnet.speciation.MulSpeciesTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */
public class MulTreeNodeSlideParser extends AbstractXMLObjectParser {
    public static final String MULTREE_NODE_REHEIGHT = "mulTreeNodeReHeight";

    public String getParserName() {
        return MULTREE_NODE_REHEIGHT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	MulSpeciesBindings species = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);
    	MulSpeciesTreeModel tree = (MulSpeciesTreeModel) xo.getChild(MulSpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
//            final double range = xo.getAttribute("range", 1.0);
//            if( range <= 0 || range > 1.0 ) {
//                throw new XMLParseException("range out of range");
//            }
        //final boolean oo = xo.getAttribute("outgroup", false);
        return new MulTreeNodeSlide(tree, species /*, range*//*, oo*/, weight);
    }

    public String getParserDescription() {
        return "Specialized Species tree operator, transform tree without breaking embedding of gene trees.";
    }

    public Class getReturnType() {
        return MulTreeNodeSlide.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
               // AttributeRule.newDoubleRule("range", true),
              //  AttributeRule.newBooleanRule("outgroup", true),

                new ElementRule(MulSpeciesBindings.class),
                new ElementRule(MulSpeciesTreeModel.class)
        };
    }

}
