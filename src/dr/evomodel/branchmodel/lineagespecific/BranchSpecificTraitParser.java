/*
 * BranchSpecificTraitParser.java
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

package dr.evomodel.branchmodel.lineagespecific;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BranchSpecificTraitParser extends AbstractXMLObjectParser {

	public static final String BRANCH_SPECIFIC_TRAIT = "branchSpecificTrait";

	@Override
	public String getParserName() {
		return BRANCH_SPECIFIC_TRAIT;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
//		CompoundParameter parameter = (CompoundParameter) xo.getChild(CompoundParameter.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		
		return new BranchSpecificTrait(treeModel, branchModel, xo.getId());
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {

		new ElementRule(BranchModel.class, false), //
				new ElementRule(TreeModel.class, false), //

		};
	}

	@Override
	public String getParserDescription() {
		return BRANCH_SPECIFIC_TRAIT;
	}

	@Override
	public Class getReturnType() {
		return BranchSpecificTrait.class;
	}

}// END: class
