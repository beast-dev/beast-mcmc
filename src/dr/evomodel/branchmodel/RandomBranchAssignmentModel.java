/*
 * RandomBranchAssignmentModel.java
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

package dr.evomodel.branchmodel;

import java.util.LinkedHashMap;
import java.util.List;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

@SuppressWarnings("serial")
public class RandomBranchAssignmentModel extends AbstractModel implements BranchModel {

	public static final boolean DEBUG = true;
	
    public static final String RANDOM_BRANCH_ASSIGNMENT_MODEL = "randomBranchAssignmentModel";
    private final TreeModel treeModel;
    private final List<SubstitutionModel> substitutionModels;
    
//    private int[] order;
    private LinkedHashMap<NodeRef, Integer> branchAssignmentMap;
    
	public RandomBranchAssignmentModel(TreeModel treeModel,
            List<SubstitutionModel> substitutionModels) {
		
		super(RANDOM_BRANCH_ASSIGNMENT_MODEL);
		
		
		this.treeModel = treeModel;
		this.substitutionModels = substitutionModels;
		
		int nodeCount = treeModel.getNodeCount();
		int nModels = substitutionModels.size();
		
		// randomly decide order, once and for all
		branchAssignmentMap = new LinkedHashMap<NodeRef, Integer>();
		for (int i = 0; i < nodeCount; i++) {

			NodeRef node = treeModel.getNode(i);
			int branchClass = Integer.MAX_VALUE; //MathUtils.nextInt(nModels);
			
			if(DEBUG) {
				
//				System.out.println(node.toString());
				
				// hack to get fixed indexing
				if(node.toString().equalsIgnoreCase("node 0, height=0.0: SimSeq1") || 
						node.toString().equalsIgnoreCase("node 1, height=0.0: SimSeq2") ||
						node.toString().equalsIgnoreCase("node 4, height=22.0")
						) {
					branchClass = 0; // 5
				} else {
					branchClass = 1; // 10
				}//END: node check
				
			} else {
				branchClass = MathUtils.nextInt(nModels);
			}//END: DEBUG check
			
			branchAssignmentMap.put(node, branchClass);
			
		}// END: nodes loop		
		
	}//END: Constructor
	
	@Override
	public Mapping getBranchModelMapping(NodeRef branch) {
		
		final int branchClass = branchAssignmentMap.get(branch);
		
        return new Mapping() {
            public int[] getOrder() {
                return new int[] { branchClass };
            }

            public double[] getWeights() {
                return new double[] { 1.0 };
            }
        };
	}

	@Override
	public List<SubstitutionModel> getSubstitutionModels() {
		return substitutionModels;
	}

	@Override
	public SubstitutionModel getRootSubstitutionModel() {
		int rootClass = branchAssignmentMap.get(treeModel.getRoot());
		return substitutionModels.get(rootClass);
	}

	@Override
	public FrequencyModel getRootFrequencyModel() {
		return getRootSubstitutionModel().getFrequencyModel();
	}

	@Override
	public boolean requiresMatrixConvolution() {
		return false;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		fireModelChanged();
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
	}

	@Override
	protected void storeState() {
	}

	@Override
	protected void restoreState() {
	}

	@Override
	protected void acceptState() {
	}


}//END: class
