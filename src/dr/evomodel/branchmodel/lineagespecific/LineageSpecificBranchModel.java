/*
 * LineageSpecificBranchModel.java
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

import java.util.*;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.codon.MG94CodonModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

/**
 * @author Filip Bielejec
 * @author Guy Baele
 * @author Marc A. Suchard
 * @version $Id$
 */
@SuppressWarnings("serial")
public class LineageSpecificBranchModel extends AbstractModel implements BranchModel, Citable {

	 public static final String LINEAGE_SPECIFIC_BRANCH_MODEL = "lineageSpecificBranchModel";
	
	private static final boolean DEBUG = false;
	
    private boolean setupMapping = true;
	
    private Map<NodeRef, Mapping> nodeMap;
    private List<SubstitutionModel> substitutionModels;
    private TreeModel treeModel;
    private FrequencyModel rootFrequencyModel;
    
    // for discrete categories
    private CountableBranchCategoryProvider categoriesProvider;
    private Parameter categoriesParameter;

	public LineageSpecificBranchModel(TreeModel treeModel, //
			FrequencyModel rootFrequencyModel, //
			final List<SubstitutionModel> substitutionModels, //
//			CountableBranchCategoryProvider categoriesProvider, //
			Parameter categoriesParameter //
	) {

		super("");

		this.treeModel = treeModel;
		this.substitutionModels = substitutionModels;
//		this.categoriesProvider = categoriesProvider;
		this.categoriesParameter = categoriesParameter;
		this.rootFrequencyModel = rootFrequencyModel;

		this.categoriesProvider = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(treeModel, categoriesParameter);
		
		this.nodeMap = new HashMap<NodeRef, Mapping>();

		for (SubstitutionModel model : this.substitutionModels) {
			addModel(model);
		}

		addModel(this.treeModel);
		addModel(this.rootFrequencyModel);
		addModel((Model) this.categoriesProvider);
		addVariable(this.categoriesParameter);

	}// END: Constructor
	
    @Override
    public Mapping getBranchModelMapping(NodeRef branch) {

        if (setupMapping) {
            setupNodeMap(branch);
        }

        return nodeMap.get(branch);
    }//END: getBranchModelMapping

	private void setupNodeMap(NodeRef branch) {

		if (branch != treeModel.getRoot()) {//TODO: neccessary?

			int branchCategory = categoriesProvider.getBranchCategory(
					treeModel, branch);
			final int uCategory = (int) categoriesParameter.getParameterValue(branchCategory);

			if (DEBUG) {
				System.out.println("branch length: " + treeModel.getBranchLength(branch) + ", " + "category:" + uCategory);
			}// END: DEBUG

			nodeMap.put(branch, new Mapping() {

				@Override
				public int[] getOrder() {
					return new int[] { uCategory };
				}// END: getOrder

				@Override
				public double[] getWeights() {
					return new double[] { 1.0 };
				}// END: getWeights

			});

		}//END: root check
		
	}// END: setupNodeMap

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
    	return substitutionModels;
    }

	@Override
	public SubstitutionModel getRootSubstitutionModel() {
		throw new RuntimeException("This method should never be called!");
	}

    @Override
    public FrequencyModel getRootFrequencyModel() {
        return rootFrequencyModel;
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
    	fireModelChanged();
    }

	@Override
	protected void storeState() {
		//
	}// END: storeState

	@Override
	protected void restoreState() {
		setupMapping = true;
	}// END: restoreState

	@Override
	protected void acceptState() {
		//
	}// END: acceptState

    public static void main(String[] args) {

        try {

            // the seed of the BEAST
            MathUtils.setSeed(666);

            // create tree
            NewickImporter importer = new NewickImporter(
                    "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
            TreeModel tree = new TreeModel(importer.importTree(null));

            // create site model
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                    "siteModel");

            // create branch rate model
            BranchRateModel branchRateModel = new DefaultBranchRateModel();

            int sequenceLength = 10;
            ArrayList<Partition> partitionsList = new ArrayList<Partition>();

            // create Frequency Model
            Parameter freqs = new Parameter.Default(new double[]{
                    0.0163936, //
                    0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
                    0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
                    0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
                    0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, //
                    0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344, 0.01639344 //
            });
            FrequencyModel freqModel = new FrequencyModel(Codons.UNIVERSAL,
                    freqs);

            // create substitution model
            Parameter alpha = new Parameter.Default(1, 10);
            Parameter beta = new Parameter.Default(1, 5);
            MG94CodonModel mg94 = new MG94CodonModel(Codons.UNIVERSAL, alpha, beta, freqModel);

            HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(mg94);

            // create partition
            Partition partition1 = new Partition(tree, //
                    substitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    0, // from
                    sequenceLength - 1, // to
                    1 // every
            );

            partitionsList.add(partition1);

            // feed to sequence simulator and generate data
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
                    partitionsList);

            Alignment alignment = simulator.simulate(false, false);

            ConvertAlignment convert = new ConvertAlignment(Nucleotides.INSTANCE,
                    GeneticCode.UNIVERSAL, alignment);


			List<SubstitutionModel> substModels = new ArrayList<SubstitutionModel>();
			for (int i = 0; i < 2; i++) {
//				alpha = new Parameter.Default(1, 10 );
//				beta = new Parameter.Default(1, 5 );
//				mg94 = new MG94CodonModel(Codons.UNIVERSAL, alpha, beta,
//						freqModel);
				substModels.add(mg94);
			}
            
			Parameter uCategories = new Parameter.Default(2, 0);
//            CountableBranchCategoryProvider provider = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(tree, uCategories);
			
            LineageSpecificBranchModel branchSpecific = new LineageSpecificBranchModel(tree, freqModel, substModels, //provider, 
            		uCategories);

            BeagleTreeLikelihood like = new BeagleTreeLikelihood(convert, //
                    tree, //
                    branchSpecific, //
                    siteRateModel, //
                    branchRateModel, //
                    null, //
                    false, //
                    PartialsRescalingScheme.DEFAULT, true);

            BeagleTreeLikelihood gold = new BeagleTreeLikelihood(convert, //
                    tree, //
                    substitutionModel, //
                    siteRateModel, //
                    branchRateModel, //
                    null, //
                    false, //
                    PartialsRescalingScheme.DEFAULT, true);
            
            System.out.println("likelihood (gold) = " + gold.getLogLikelihood());
            System.out.println("likelihood = " + like.getLogLikelihood());
            
        } catch (Exception e) {
            e.printStackTrace();
        }

    }// END: main

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.MOLECULAR_CLOCK;
	}

	@Override
	public String getDescription() {
		return "Lineage Specific Branch model";
	}

	public List<Citation> getCitations() {
		return Collections.singletonList(
				new Citation(new Author[] { new Author("F", "Bielejec"),
				new Author("P", "Lemey"), new Author("G", "Baele"), new Author("A", "Rambaut"),
				new Author("MA", "Suchard") }, Citation.Status.IN_PREPARATION));
	}// END: getCitations

}// END: class
