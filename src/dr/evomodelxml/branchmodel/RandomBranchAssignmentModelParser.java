/*
 * RandomBranchAssignmentModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchmodel;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import dr.evomodel.branchmodel.RandomBranchAssignmentModel;
import dr.evomodel.branchmodel.RandomBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class RandomBranchAssignmentModelParser extends AbstractXMLObjectParser {

	 public static final String MODELS = "models";
	
	@Override
	public String getParserName() {
		return RandomBranchAssignmentModel.RANDOM_BRANCH_ASSIGNMENT_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("\nUsing random assignment branch model.");
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

	      XMLObject cxo = xo.getChild(MODELS);
	      List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
	      for (int i = 0; i < cxo.getChildCount(); i++) {

				SubstitutionModel substModel = (SubstitutionModel) cxo.getChild(i);
				substitutionModels.add(substModel);
				
			}//END: models loop
		
		return new RandomBranchAssignmentModel(treeModel, substitutionModels);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		
		return new XMLSyntaxRule[]{
				
				new ElementRule(TreeModel.class, false), //
                new ElementRule(MODELS,
                        new XMLSyntaxRule[] {
                                new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE),
                        }
                )             
                
		};
	}//END: XMLSyntaxRule

	@Override
	public String getParserDescription() {
		return "This element provides a branch model which randomly assigns " +
				"substitution models to branches on the tree by sampling " +
				"with replacement from the provided list of substitution models. ";
	}

	@Override
	public Class getReturnType() {
		return RandomBranchModel.class;
	}


}//END: class
