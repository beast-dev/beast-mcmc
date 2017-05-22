/*
 * BranchSpecificTrait.java
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

import dr.evolution.tree.*;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BranchSpecificTrait implements TreeTraitProvider {

	private Helper helper;
	private TreeModel treeModel;
	
	public BranchSpecificTrait(
		TreeModel treeModel,
		Parameter uCategories,
		CompoundParameter uniqueParameters, 
		final String parameterName) {
		
		this.treeModel = treeModel;
		helper = new Helper();
		
		TreeTrait<Double> uTrait = new TreeTrait.D() {

			@Override
			public String getTraitName() {
				return parameterName;
			}

			@Override
			public dr.evolution.tree.TreeTrait.Intent getIntent() {
				return Intent.BRANCH;
			}

			@Override
			public Double getTrait(Tree tree, NodeRef node) {

				double value = 0.0;
				
				
				
				
				
				
				
				
				
				
				
				return null;
			}
			
			
			
		};
		
		
	}//END: Constructor
	
	public BranchSpecificTrait(
			TreeModel treeModel,
			final BranchModel branchModel,
//			, final CompoundParameter parameter 
			final String parameterName
			) {
		
		this.treeModel = treeModel;
		helper = new Helper();
		
		//TODO: this could annotate with all Variables in Substitution model
		TreeTrait<Double> uTrait = new TreeTrait.D() {

			@Override
			public String getTraitName() {
				return parameterName;//parameter.getId();
			}

			@Override
			public dr.evolution.tree.TreeTrait.Intent getIntent() {
				return Intent.BRANCH;
			}

			@Override
			public Double getTrait(Tree tree, NodeRef branch) {

				double value = 0.0;
				
				int[] uCats = branchModel.getBranchModelMapping(branch).getOrder();
				int category = uCats[0];

				SubstitutionModel substmodel = branchModel.getSubstitutionModels().get(category);
                value = (Double) substmodel.getVariable(0).getValue(0);		
				
				
				return value;

			}//END: getTrait
		};
		
		helper.addTrait(uTrait);
		
	}//END: Constructor

    public TreeTrait[] getTreeTraits() {
        return helper.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return helper.getTreeTrait(key);
    }

	public String toString() {

		String annotatedTree = TreeUtils.newick(treeModel,
				new TreeTraitProvider[] { this });

		return annotatedTree;

	}// END: toString
    
}//END: class
