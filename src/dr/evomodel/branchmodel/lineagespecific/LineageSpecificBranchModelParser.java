/*
 * LineageSpecificBranchModelParser.java
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

import java.util.ArrayList;
import java.util.List;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Parameter;
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
public class LineageSpecificBranchModelParser extends AbstractXMLObjectParser {

	 public static final String MODELS = "models";
	 public static final String CATEGORIES = "categories";
	
	@Override
	public String getParserName() {
		return LineageSpecificBranchModel.LINEAGE_SPECIFIC_BRANCH_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		FrequencyModel rootFrequencyModel = (FrequencyModel) xo.getChild(FrequencyModel.class);
		
	      List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
	      
	      XMLObject cxo = xo.getChild(MODELS);
			for (int i = 0; i < cxo.getChildCount(); i++) {

				SubstitutionModel substModel = (SubstitutionModel) cxo.getChild(i);
				substitutionModels.add(substModel);
				
			}//END: models loop
		
			// TODO: check if categories numbering starts from zero
			Parameter categories = (Parameter) xo.getElementFirstChild(CATEGORIES); 
			
//			CountableBranchCategoryProvider.CladeBranchCategoryModel provider = new CountableBranchCategoryProvider.CladeBranchCategoryModel(treeModel, categories);// false);
			
		return new LineageSpecificBranchModel(treeModel, rootFrequencyModel, substitutionModels, //provider, 
				categories);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				
				new ElementRule(TreeModel.class, false), //
				new ElementRule(FrequencyModel.class, false), //
                new ElementRule(MODELS,
                        new XMLSyntaxRule[] {
                                new ElementRule(AbstractModel.class, 1, Integer.MAX_VALUE),
                        }
                ), //
                new ElementRule(CATEGORIES,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class, 1, 1),
                        }
                ) //              
                
                
		};
	}

	@Override
	public String getParserDescription() {
		return "This element provides a branch model which has branches assigned to specific substitution models." +
				"These assignments can then be changed in course of MCMC.";
	}

	@Override
	public Class getReturnType() {
		return LineageSpecificBranchModel.class;
	}

}//END: class
