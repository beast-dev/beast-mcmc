/*
 * RandomBranchModel.java
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
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math.random.MersenneTwister;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */

@SuppressWarnings("serial")
public class RandomBranchModel extends AbstractModel implements BranchModel {

	public static final String RANDOM_BRANCH_MODEL = "randomBranchModel";
	private final TreeModel treeModel;
	private GY94CodonModel baseSubstitutionModel;
	private LinkedList<SubstitutionModel> substitutionModels;
	private LinkedHashMap<NodeRef, Integer> branchAssignmentMap;

	// TODO: parse distribution model for e_i
	// TODO: parse substitution model (hardcoded now)
	// TODO: parse parameter that follows trend (hardcoded now)
	private double rate;
	
	private static MersenneTwister random;
	
	public RandomBranchModel(TreeModel treeModel, //
			GY94CodonModel baseSubstitutionModel, //
			double rate, //
			boolean hasSeed,
			long seed
			) {

		super(RANDOM_BRANCH_MODEL);

		this.treeModel = treeModel;
		this.baseSubstitutionModel = baseSubstitutionModel;

		this.rate = rate;

		if( hasSeed) {
			
			// use fixed seed for e_i
			random = new MersenneTwister(
					seed);
		} else {

			//use BEAST seed
			random = new MersenneTwister( MathUtils.nextLong() );
		
		}//END: seed check
		
		
		setup();

	}// END: Constructor

	private void setup() {

		DataType dataType = baseSubstitutionModel.getDataType();
		FrequencyModel freqModel = baseSubstitutionModel.getFrequencyModel();
		Parameter kappaParameter = new Parameter.Default("kappa", 1,
				baseSubstitutionModel.getKappa());

		substitutionModels = new LinkedList<SubstitutionModel>();
		branchAssignmentMap = new LinkedHashMap<NodeRef, Integer>();

		int branchClass = 0;
		for (NodeRef node : treeModel.getNodes()) {
			if (!treeModel.isRoot(node)) {

				double nodeHeight = treeModel.getNodeHeight(node);
				double parentHeight = treeModel.getNodeHeight(treeModel
						.getParent(node));

				double time = 0.5 * (parentHeight + nodeHeight);
				
				double baseOmega = baseSubstitutionModel.getOmega();

				double fixed = baseOmega * time;
				
				double epsilon = (Math.log(1-random.nextDouble()) / (-rate) );  //Math.exp((random.nextGaussian() * stdev + mean));
				
				double value = fixed + epsilon;
				
				Parameter omegaParameter = new Parameter.Default("omega", 1,
						value);

				GY94CodonModel gy94 = new GY94CodonModel((Codons) dataType,
						omegaParameter, kappaParameter, freqModel);

				substitutionModels.add(gy94);
				branchAssignmentMap.put(node, branchClass);
				branchClass++;
			}//END: root check
		}// END: nodes loop

	}// END: setup

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
		// int rootClass = branchAssignmentMap.get(treeModel.getRoot());
		// return substitutionModels.get(rootClass);
		throw new RuntimeException("Not implemented!");
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

}// END: class
