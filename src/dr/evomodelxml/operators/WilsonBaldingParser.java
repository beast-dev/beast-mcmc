/*
 * WilsonBaldingParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeOperator;
import dr.evomodel.treelikelihood.thorneytreelikelihood.UniformSubtreePruneRegraft;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class WilsonBaldingParser extends AbstractXMLObjectParser {

    public static final String WILSON_BALDING = "wilsonBalding";
    public static final String DEMOGRAPHIC_MODEL = "demographicModel";

    public String getParserName() {
        return WILSON_BALDING;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        WilsonBalding op = new WilsonBalding(treeModel,weight);
        if(treeModel instanceof ConstrainedTreeModel){
            return new ConstrainedTreeOperator((ConstrainedTreeModel) treeModel, weight, op);
        }else{
            return op;
        }
        
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };

    public String getParserDescription() {
        return "An operator which performs the Wilson-Balding move on a tree";
    }

    public Class getReturnType() {
        return WilsonBalding.class;
    }
}
