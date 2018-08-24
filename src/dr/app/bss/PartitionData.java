/*
 * PartitionData.java
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

package dr.app.bss;

import java.io.Serializable;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.codon.MG94HKYCodonModel;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.substmodel.aminoacid.*;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.TaxaParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.InverseGaussianDistributionModel;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inferencexml.distribution.DistributionModelParser;
import dr.inferencexml.distribution.InverseGaussianDistributionModelParser;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionData implements Serializable {

	public PartitionData() {
	}// END: Constructor

	public int from = 1;
	public int to = 10;
	public int every = 1;

	public int createPartitionSiteCount() {
		return ((to - from) / every) + 1;
	}

	public void resetIdrefs() {
		resetClockModelIdref();
		resetFrequencyModelIdref();
		resetSiteRateModelIdref();
		resetSubstitutionModelIdref();
		resetTreeModelIdref();
		resetDemographicModelIdref();
		resetTaxaIdref();
	}

	// ///////////////////////
	// ---TREE ANNOTATING---//
	// ///////////////////////
	
//	private LinkedHashMap<NodeRef, int[]> sequenceMap;
//	
//	public void setSequenceMap(LinkedHashMap<NodeRef, int[]> sequenceMap) {
//		this.sequenceMap = sequenceMap;
//	}
//	
//	public LinkedHashMap<NodeRef, int[]> getSequenceMap() {
//		return sequenceMap;
//	}
	
	// /////////////////////////
	// ---DEMOGRAPHIC MODEL---//
	// /////////////////////////

	//TODO: LogisticGrowth.getInverseIntensity
	public static final int lastImplementedIndex = 3;
	
	public int demographicModelIndex = 0;
	
	public String demographicModelIdref = Utils.DEMOGRAPHIC_MODEL;

	public void resetDemographicModelIdref() {
		this.demographicModelIdref = Utils.DEMOGRAPHIC_MODEL;
	}
	
    public static String[] demographicModels = {
    	"No Model (user-specified tree)",
    	"Constant Population",
        "Exponential Growth (Growth Rate)",
        "Exponential Growth (Doubling Time)",
//        "Logistic Growth (Growth Rate)",
//        "Logistic Growth (Doubling Time)",
//        "Expansion (Growth Rate)",
//        "Expansion (Doubling Time)",
        };
	
	public static String[] demographicParameterNames = new String[] {
	        "Population Size", // Constant Population
			
	        "Population Size", // Exponential Growth (Growth Rate)
			"Growth Rate", // Exponential Growth (Growth Rate)
			
			"Population Size", // Exponential Growth (Doubling Time)
			"Doubling Time", // Exponential Growth (Doubling Time)
			
//			"Population Size", // Logistic Growth (Growth Rate)
//			"Growth Rate", // Logistic Growth (Growth Rate)
//			"Logistic Shape (Half-life)", // Logistic Growth (Growth Rate)
			
//			"Population Size", // Logistic Growth (Doubling Time)
//			"Doubling Time", // Logistic Growth (Doubling Time)
//			"Logistic Shape (Half-life)", // Logistic Growth (Doubling Time)
			
//			"Population Size", // Expansion (Growth Rate)
//			"Ancestral Proportion", // Expansion (Growth Rate)
//			"Growth Rate", // Expansion (Growth Rate)
			
//			"Population Size", // Expansion (Doubling Time)
//			"Ancestral Proportion", // Expansion (Doubling Time)
//			"Doubling Time", // Expansion (Doubling Time)
	};	
	
	public static int[][] demographicParameterIndices = {
		    {  }, // No model
			{ 0 }, // Constant Population
			{ 1, 2 }, // Exponential Growth (Growth Rate)
			{ 3, 4 }, // Exponential Growth (Doubling Time)
//			{ 5, 6, 7 }, // Logistic Growth (Growth Rate)
//			{ 8, 9, 10 }, // Logistic Growth (Doubling Time)
//			{ 11, 12, 13 }, // Expansion (Growth Rate)
//			{ 14, 15, 16 } // Expansion (Doubling Time)
	};
	
	
	public double[] demographicParameterValues = new double[] {
			
			/*Constant Population*/
			1000.0, // Population Size
			
			/*Exponential Growth (Growth Rate)*/
			1000.0, // Population Size
			0.5, // Growth Rate
			
			/*Exponential Growth (Doubling Time)*/
			1000.0, // Population Size
			10.0, // Doubling Time
			
//			/*Logistic Growth (Growth Rate)*/
//			1000.0, // Population Size
//			0.5, // Growth Rate
//			50.0, // Logistic Shape (Half-life)
			
//			1000.0, // Population Size
//			10.0, // Doubling Time
//			50.0, // Logistic Shape (Half-life)
//			1000.0, // Population Size
//			0.1, // Ancestral Proportion
//			0.5, // Growth Rate
//			1000.0, // Population Size
//			0.1, // Ancestral Proportion
//			10.0 // Doubling Time
	};

	public DemographicFunction createDemographicFunction() {

		DemographicFunction demographicFunction = null;

		if (this.demographicModelIndex == 0) { // No model

			// do nothing
			
		} else if (this.demographicModelIndex == 1) {// Constant Population

			demographicFunction = new ConstantPopulation(Units.Type.YEARS);
			((ConstantPopulation)demographicFunction).setN0(demographicParameterValues[0]);
			
		} else if (this.demographicModelIndex == 2) {// Exponential Growth (Growth Rate)

			demographicFunction = new ExponentialGrowth(Units.Type.YEARS);
            ((ExponentialGrowth) demographicFunction).setN0(demographicParameterValues[1]);
            ((ExponentialGrowth) demographicFunction).setGrowthRate(demographicParameterValues[2]);
			
		} else if (this.demographicModelIndex == 3) {// Exponential Growth (Doubling Time)
			
			demographicFunction = new ExponentialGrowth(Units.Type.YEARS);
            ((ExponentialGrowth) demographicFunction).setN0(demographicParameterValues[3]);
            ((ExponentialGrowth) demographicFunction).setDoublingTime(demographicParameterValues[4]);
			
//		} else if (this.demographicModelIndex == 4) {// Logistic Growth (Growth Rate)
//			
//			demographicFunction = new LogisticGrowth(Units.Type.YEARS);
//            ((LogisticGrowth) demographicFunction).setN0(demographicParameterValues[5]);
//            ((LogisticGrowth) demographicFunction).setGrowthRate(demographicParameterValues[6]);
//            ((LogisticGrowth) demographicFunction).setTime50(demographicParameterValues[7]);
//			
//		} else if (this.demographicModelIndex == 5) {// Logistic Growth (Doubling Time)
//			
//			demographicFunction = new LogisticGrowth(Units.Type.YEARS);
//            ((LogisticGrowth) demographicFunction).setN0(demographicParameterValues[8]);
//            ((LogisticGrowth) demographicFunction).setDoublingTime(demographicParameterValues[9]);
//            ((LogisticGrowth) demographicFunction).setTime50(demographicParameterValues[10]);
//			
//		} else if (this.demographicModelIndex == 6) {// Expansion (Growth Rate)
//			
//			demographicFunction = new Expansion(Units.Type.YEARS);
//            ((Expansion) demographicFunction).setN0(demographicParameterValues[11]);
//            ((Expansion) demographicFunction).setProportion(demographicParameterValues[12]);
//            ((Expansion) demographicFunction).setGrowthRate(demographicParameterValues[13]);
//			
//		} else if (this.demographicModelIndex == 7) {// Expansion (Doubling Time)
//			
//			demographicFunction = new Expansion(Units.Type.YEARS);
//            ((Expansion) demographicFunction).setN0(demographicParameterValues[14]);
//            ((Expansion) demographicFunction).setProportion(demographicParameterValues[15]);
//            ((Expansion) demographicFunction).setDoublingTime(demographicParameterValues[16]);
			
		} else {

			System.out.println("Not yet implemented");
			
		}

		return demographicFunction;
	}// END: createDemographicFunction
	
	// ////////////
	// ---TAXA---//
	// ////////////

	public TreesTableRecord record = null;
	
	public String taxaIdref = TaxaParser.TAXA;
	
	public void resetTaxaIdref() {
		this.taxaIdref = TaxaParser.TAXA;
	}
	
	// //////////////////
	// ---TREE MODEL---//
	// //////////////////

//	public Tree tree = null;
	public String treeModelIdref = TreeModel.TREE_MODEL;

	public void resetTreeModelIdref() {
	this.treeModelIdref = TreeModel.TREE_MODEL;
	}
	
	public TreeModel createTreeModel() {
		
		TreeModel treeModel = null;
		if (this.demographicModelIndex == 0 && this.record.isTreeSet()) {
			
			treeModel = new TreeModel(this.record.getTree());
			
		} else if( (this.demographicModelIndex > 0 && this.demographicModelIndex <= lastImplementedIndex) && this.record.isTreeSet()) {
			
			Taxa taxa = new Taxa(this.record.getTree().asList()); 
			CoalescentSimulator topologySimulator = new CoalescentSimulator();
			treeModel = new TreeModel(topologySimulator.simulateTree(taxa, createDemographicFunction()));			
			
		} else if((this.demographicModelIndex > 0 && this.demographicModelIndex <= lastImplementedIndex) && this.record.isTaxaSet()) {
			
			Taxa taxa = this.record.getTaxa();
			CoalescentSimulator topologySimulator = new CoalescentSimulator();
			treeModel = new TreeModel(topologySimulator.simulateTree(taxa, createDemographicFunction()));
		
//			} else if (this.demographicModelIndex == 0 && this.record.taxaSet) { 
//			throw new RuntimeException("Data and demographic model incompatible for partition ");	
			
		} else {
			
			throw new RuntimeException("Data and demographic model incompatible.");
			
		}// END: demo model check
		
		return treeModel;
	}//END: createTreeModel

	// /////////////////
	// ---DATA TYPE---//
	// /////////////////

	public int dataTypeIndex = 0;

	public static String[] dataTypes = { "Nucleotide", //
			"Codon", //
			"Amino acid" //
	};

	public DataType createDataType() {

		DataType dataType = null;

		if (this.dataTypeIndex == 0) { // Nucleotide

			dataType = Nucleotides.INSTANCE;

		} else if (this.dataTypeIndex == 1) { // Codon

			dataType = Codons.UNIVERSAL;

		} else if (this.dataTypeIndex == 2) { // AminoAcid

			dataType = AminoAcids.INSTANCE;
			
		} else {

			System.out.println("Not yet implemented");

		}

		return dataType;
	}// END: createDataType

	// ///////////////////////////
	// ---SUBSTITUTION MODELS---//
	// ///////////////////////////

	public int substitutionModelIndex = 0;

	public String substitutionModelIdref = Utils.SUBSTITUTION_MODEL;

	public void resetSubstitutionModelIdref() {
		this.substitutionModelIdref = Utils.SUBSTITUTION_MODEL;
	}
	
	public static String[] substitutionModels = { "HKY", //
			"GTR", //
			"TN93", //
			"GY94CodonModel", //
			"MG94CodonModel",
            "Blosum62", //	
			"CPREV", //
			"Dayhoff", //
			"FLU", //
			"JTT", //
			"LG", //
			"MTREV", //
			"WAG" //
	};

	public static int[] substitutionCompatibleDataTypes = { 0, // HKY
			0, // GTR
			0, // TN93
			1, // GY94CodonModel
			1, // MG94CodonModel
			2, // Blosum62
			2, // CPREV
			2, // Dayhoff
			2, // FLU
			2, // JTT
			2, // LG
			2, // MTREV
			2 // WAG
	};
	
	public static String[] substitutionParameterNames = new String[] {
			"Kappa value", // HKY
			"AC", // GTR
			"AG", // GTR
			"AT", // GTR
			"CG", // GTR
			"CT", // GTR
			"GT", // GTR
			"Kappa 1 (A-G)", // TN93
			"Kappa 2 (C-T)", // TN93
			"Omega value", // GY94CodonModel
			"Kappa value", // GY94CodonModel
			"Alpha value", // MG94CodonModel
			"Beta value", // MG94CodonModel
			"Kappa value" // MG94CodonModel
			
	};

	public static int[][] substitutionParameterIndices = { { 0 }, // HKY
			{ 1, 2, 3, 4, 5, 6 }, // GTR
			{ 7, 8 }, // TN93
			{ 9, 10 }, // GY94CodonModel
			{11, 12, 13}, // MG94CodonModel
			{}, // Blosum62
			{}, // CPREV
			{}, // Dayhoff
			{}, // FLU
			{}, // JTT
			{}, // LG
			{}, // MTREV
			{} // WAG
	};

	public double[] substitutionParameterValues = new double[] { 1.0, // Kappa-value
			1.0, // AC
			1.0, // AG
			1.0, // AT
			1.0, // CG
			1.0, // CT
			1.0, // GT
			1.0, // Kappa 1
			1.0, // Kappa 2
			0.1, // Omega value
			1.0, // Kappa value
			10.0, // Alpha value
			1.0, // Beta value
			1.0 // kappa value
	};

	public BranchModel createBranchModel() {

		BranchModel branchModel = null;

		if (this.substitutionModelIndex == 0) { // HKY

			Parameter kappa = new Parameter.Default(1, substitutionParameterValues[0]);

			FrequencyModel frequencyModel = this.createFrequencyModel();

			HKY hky = new HKY(kappa, frequencyModel);

			branchModel = new HomogeneousBranchModel(hky);

		} else if (this.substitutionModelIndex == 1) { // GTR

			Parameter ac = new Parameter.Default(1,
					substitutionParameterValues[1]);
			Parameter ag = new Parameter.Default(1,
					substitutionParameterValues[2]);
			Parameter at = new Parameter.Default(1,
					substitutionParameterValues[3]);
			Parameter cg = new Parameter.Default(1,
					substitutionParameterValues[4]);
			Parameter ct = new Parameter.Default(1,
					substitutionParameterValues[5]);
			Parameter gt = new Parameter.Default(1,
					substitutionParameterValues[6]);

			FrequencyModel frequencyModel = this.createFrequencyModel();

			GTR gtr = new GTR(ac, ag, at, cg, ct, gt, frequencyModel);

			branchModel = new HomogeneousBranchModel(gtr);

		} else if (this.substitutionModelIndex == 2) { // TN93

			Parameter kappa1 = new Parameter.Default(1,
					substitutionParameterValues[7]);
			Parameter kappa2 = new Parameter.Default(1,
					substitutionParameterValues[8]);

			FrequencyModel frequencyModel = this.createFrequencyModel();

			TN93 tn93 = new TN93(kappa1, kappa2, frequencyModel);

			branchModel = new HomogeneousBranchModel(tn93);

		} else if (this.substitutionModelIndex == 3) { // Yang Codon Model

			FrequencyModel frequencyModel = this.createFrequencyModel();

			Parameter kappa = new Parameter.Default(1,
					substitutionParameterValues[9]);
			Parameter omega = new Parameter.Default(1,
					substitutionParameterValues[10]);

			GY94CodonModel yangCodonModel = new GY94CodonModel(
					Codons.UNIVERSAL, omega, kappa, frequencyModel);

			branchModel = new HomogeneousBranchModel(yangCodonModel);

			
		} else if(this.substitutionModelIndex == 4) { // MG94CodonModel
			
			
			FrequencyModel frequencyModel = this.createFrequencyModel();

			Parameter alpha = new Parameter.Default(1, substitutionParameterValues[11]);
			Parameter beta = new Parameter.Default(1, substitutionParameterValues[12]);
			Parameter kappa = new Parameter.Default(1, substitutionParameterValues[13]);
			
			MG94HKYCodonModel mg94 = new MG94HKYCodonModel(Codons.UNIVERSAL, alpha, beta, kappa, frequencyModel);

			branchModel = new HomogeneousBranchModel(mg94);
			
		} else if (this.substitutionModelIndex == 5) { // Blosum62
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = Blosum62.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);
			
        } else if (this.substitutionModelIndex == 6) { // CPREV
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = CPREV.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);
			
        } else if (this.substitutionModelIndex == 7) { // Dayhoff
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = Dayhoff.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);
		
        } else if (this.substitutionModelIndex == 8) { // JTT
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = JTT.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);
		
        } else if (this.substitutionModelIndex == 9) { // LG
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = LG.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);
			
        } else if (this.substitutionModelIndex == 10) { // MTREV
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = MTREV.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);	
			
        } else if (this.substitutionModelIndex == 11) { // WAG
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			EmpiricalRateMatrix rateMatrix = WAG.INSTANCE;

			EmpiricalAminoAcidModel empiricalAminoAcidModel = new EmpiricalAminoAcidModel(
					rateMatrix, frequencyModel);

			branchModel = new HomogeneousBranchModel(
					empiricalAminoAcidModel);	
			
		} else {

			System.out.println("Not yet implemented");

		}

		return branchModel;
	}// END: createBranchSubstitutionModel

	// ////////////////////////
	// ---FREQUENCY MODELS---//
	// ////////////////////////

	public String frequencyModelIdref = Utils.FREQUENCY_MODEL;

	public void resetFrequencyModelIdref() {
		this.frequencyModelIdref = Utils.FREQUENCY_MODEL;
	}
	
	public int frequencyModelIndex = 0;

	public static String[] frequencyModels = { "Nucleotide frequencies", //
			"Codon frequencies", //
			"Amino acid frequencies" 
			};

	public static int[] frequencyCompatibleDataTypes = { 0, // Nucleotide
			1, // Codon
			2 // Amino acid
	};

	public static String[] frequencyParameterNames = new String[] {
			"A frequency", // Nucleotide frequencies
			"C frequency", // Nucleotide frequencies
			"G frequency", // Nucleotide frequencies
			"T frequency", // Nucleotide frequencies
            "AAA frequency", // Codon frequencies
            "AAC frequency", //
            "AAG frequency", //
            "AAT frequency", //
            "ACA frequency", //
            "ACC frequency", //
            "ACG frequency", //
            "ACT frequency", //
            "AGA frequency", //
            "AGC frequency", //
            "AGG frequency", //
            "AGT frequency", //
            "ATA frequency", //
            "ATC frequency", //
            "ATG frequency", //
            "ATT frequency", //
            "CAA frequency", //
            "CAC frequency", //
            "CAG frequency", //
            "CAT frequency", //
            "CCA frequency", //
            "CCC frequency", //
            "CCG frequency", //
            "CCT frequency", //
            "CGA frequency", //
            "CGC frequency", //
            "CGG frequency", //
            "CGT frequency", //
            "CTA frequency", //
            "CTC frequency", //
            "CTG frequency", //
            "CTT frequency", //
            "GAA frequency", //
            "GAC frequency", //
            "GAG frequency", //
            "GAT frequency", //
            "GCA frequency", //
            "GCC frequency", //
            "GCG frequency", //
            "GCT frequency", //
            "GGA frequency", //
            "GGC frequency", //
            "GGG frequency", //
            "GGT frequency", //
            "GTA frequency", //
            "GTC frequency", //
            "GTG frequency", //
            "GTT frequency", //
            "TAC frequency", //
            "TAT frequency", //
            "TCA frequency", //
            "TCC frequency", //
            "TCG frequency", //
            "TCT frequency", //
            "TGC frequency", //
            "TGG frequency", //
            "TGT frequency", //
            "TTA frequency", //
            "TTC frequency", //
            "TTG frequency", //
            "TTT frequency", // Codon frequencies
            "amino acid frequency 1", //
            "amino acid frequency 2", //
            "amino acid frequency 3", //
            "amino acid frequency 4", //
            "amino acid frequency 5", //
            "amino acid frequency 6", //
            "amino acid frequency 7", //
            "amino acid frequency 8", //
            "amino acid frequency 9", //
            "amino acid frequency 10", //
            "amino acid frequency 11", //
            "amino acid frequency 12", //
            "amino acid frequency 13", //
            "amino acid frequency 14", //
            "amino acid frequency 15", //
            "amino acid frequency 16", //
            "amino acid frequency 17", //
            "amino acid frequency 18", //
            "amino acid frequency 19", //
            "amino acid frequency 20" //
			
	};

	public int[][] frequencyParameterIndices = { { 0, 1, 2, 3 }, // NucleotideFrequencies
			{ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, //
					14, 15, 16, 17, 18, 19, 20, 21, 22, 23, //
					24, 25, 26, 27, 28, 29, 30, 31, 32, 33, //
					34, 35, 36, 37, 38, 39, 40, 41, 42, 43, //
					44, 45, 46, 47, 48, 49, 50, 51, 52, 53, //
					54, 55, 56, 57, 58, 59, 60, 61, 62, 63, //
					64 }, // CodonFrequencies
			{ 65, 66, 67, 68, 69, //
					70, 71, 72, 73, 74, //
					75, 76, 77, 78, 79, //
					80, 81, 82, 83, 84 //
			} // AminoAcidFrequencies
	};

	public double[] frequencyParameterValues = new double[] { 0.25, // A frequency
			0.25, // C frequency
			0.25, // G frequency
			0.25, // T frequency
			0.0163936, // AAA frequency
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, //
			0.01639344, // TTT frequency
			0.05, // aminoacidfrequency1
			0.05, // aminoacidfrequency2
			0.05, // aminoacidfrequency3
			0.05, // aminoacidfrequency4
			0.05, // aminoacidfrequency5
			0.05, // aminoacidfrequency6
			0.05, // aminoacidfrequency7
			0.05, // aminoacidfrequency8
			0.05, // aminoacidfrequency9
			0.05, // aminoacidfrequency10
			0.05, // aminoacidfrequency11
			0.05, // aminoacidfrequency12
			0.05, // aminoacidfrequency13
			0.05, // aminoacidfrequency14
			0.05, // aminoacidfrequency15
			0.05, // aminoacidfrequency16
			0.05, // aminoacidfrequency17
			0.05, // aminoacidfrequency18
			0.05, // aminoacidfrequency19
			0.05 // aminoacidfrequency20
	};

	public FrequencyModel createFrequencyModel() {

		FrequencyModel frequencyModel = null;

		if (this.frequencyModelIndex == 0) { // Nucleotidefrequencies

			Parameter freqs = new Parameter.Default(new double[] {
					frequencyParameterValues[0], frequencyParameterValues[1],
					frequencyParameterValues[2], frequencyParameterValues[3] });

			DataType dataType = this.createDataType();
			
			frequencyModel = new FrequencyModel(dataType, freqs);

		} else if (this.frequencyModelIndex == 1) {

			Parameter freqs = new Parameter.Default(new double[] {
					frequencyParameterValues[4], frequencyParameterValues[5],
					frequencyParameterValues[6], frequencyParameterValues[7],
					frequencyParameterValues[8], frequencyParameterValues[9],
					frequencyParameterValues[10], frequencyParameterValues[11],
					frequencyParameterValues[12], frequencyParameterValues[13],
					frequencyParameterValues[14], frequencyParameterValues[15],
					frequencyParameterValues[16], frequencyParameterValues[17],
					frequencyParameterValues[18], frequencyParameterValues[19],
					frequencyParameterValues[20], frequencyParameterValues[21],
					frequencyParameterValues[22], frequencyParameterValues[23],
					frequencyParameterValues[24], frequencyParameterValues[25],
					frequencyParameterValues[26], frequencyParameterValues[27],
					frequencyParameterValues[28], frequencyParameterValues[29],
					frequencyParameterValues[30], frequencyParameterValues[31],
					frequencyParameterValues[32], frequencyParameterValues[33],
					frequencyParameterValues[34], frequencyParameterValues[35],
					frequencyParameterValues[36], frequencyParameterValues[37],
					frequencyParameterValues[38], frequencyParameterValues[39],
					frequencyParameterValues[40], frequencyParameterValues[41],
					frequencyParameterValues[42], frequencyParameterValues[43],
					frequencyParameterValues[44], frequencyParameterValues[45],
					frequencyParameterValues[46], frequencyParameterValues[47],
					frequencyParameterValues[48], frequencyParameterValues[49],
					frequencyParameterValues[50], frequencyParameterValues[51],
					frequencyParameterValues[52], frequencyParameterValues[53],
					frequencyParameterValues[54], frequencyParameterValues[55],
					frequencyParameterValues[56], frequencyParameterValues[57],
					frequencyParameterValues[58], frequencyParameterValues[59],
					frequencyParameterValues[60], frequencyParameterValues[61],
					frequencyParameterValues[62], frequencyParameterValues[63],
					frequencyParameterValues[64] });

			DataType dataType = this.createDataType();
			
			frequencyModel = new FrequencyModel(dataType, freqs);

		} else if (this.frequencyModelIndex == 2) {

			Parameter freqs = new Parameter.Default(new double[] {
					frequencyParameterValues[65], frequencyParameterValues[66], frequencyParameterValues[67], frequencyParameterValues[68],
					frequencyParameterValues[69], frequencyParameterValues[70], frequencyParameterValues[71], frequencyParameterValues[72],
					frequencyParameterValues[73], frequencyParameterValues[74], frequencyParameterValues[75], frequencyParameterValues[76],
					frequencyParameterValues[77], frequencyParameterValues[78], frequencyParameterValues[79], frequencyParameterValues[80],
					frequencyParameterValues[81], frequencyParameterValues[82], frequencyParameterValues[83], frequencyParameterValues[84]
							});

			DataType dataType = this.createDataType();
			
			frequencyModel = new FrequencyModel(dataType, freqs);

		} else {

			System.out.println("Not yet implemented");

		}

		return frequencyModel;
	}// END: createFrequencyModel
	
	// ////////////////////
	// ---CLOCK MODELS---//
	// ////////////////////
	
	public int clockModelIndex = 0;
    public int LRC_INDEX = 1;
	
	public String clockModelIdref = BranchRateModel.BRANCH_RATES;

	public void resetClockModelIdref() {
		this.clockModelIdref = BranchRateModel.BRANCH_RATES;
	}
	
	public static String[] clockModels = { "Strict Clock", // 0
			"Lognormal relaxed clock (Uncorrelated)", // 1
			"Exponential relaxed clock (Uncorrelated)", // 2
			"Inverse Gaussian relaxed clock" // 3
	};

	public static String[] clockParameterNames = new String[] { "clock.rate", // StrictClock
			"ucld.mean", // Lognormal relaxed clock
			"ucld.stdev", // Lognormal relaxed clock
			"ucld.offset", // Lognormal relaxed clock
			"uced.mean", // Exponential relaxed clock
			"uced.offset", // Exponential relaxed clock
			"ig.mean", // Inverse Gaussian
			"ig.stdev", // Inverse Gaussian
			"ig.offset" // Inverse Gaussian
	};

	public static int[][] clockParameterIndices = { { 0 }, // StrictClock
			{ 1, 2, 3 }, // Lognormal relaxed clock
			{ 4, 5 }, // Exponential relaxed clock
			{ 6, 7, 8 } // Inverse Gaussian
	};

	public double[] clockParameterValues = new double[] { 1.0, // clockrate
			0.001, // ucld.mean
			1.0, // ucld.stdev
			0.0, // ucld.offset
			1.0, // uced.mean
            0.0, // uced.offset
			0.0, // ig.mean
			1.0, // ig.stdev
			0.0 // ig.offset
	};

	public boolean lrcParametersInRealSpace = true;
	
	public BranchRateModel createClockRateModel() {

		BranchRateModel branchRateModel = null;

		if (this.clockModelIndex == 0) { // Strict Clock

			Parameter rateParameter = new Parameter.Default(1, clockParameterValues[0]);
			
			branchRateModel = new StrictClockBranchRates(rateParameter);

		} else if (this.clockModelIndex == LRC_INDEX) {// Lognormal relaxed clock

			double numberOfBranches = 2 * (createTreeModel().getTaxonCount() - 1);
			Parameter rateCategoryParameter = new Parameter.Default(numberOfBranches);
			
			Parameter mean = new Parameter.Default(LogNormalDistributionModelParser.MEAN, 1, clockParameterValues[1]);
			Parameter stdev = new Parameter.Default(LogNormalDistributionModelParser.STDEV, 1, clockParameterValues[2]);
			//TODO: choose between log scale / real scale
	        ParametricDistributionModel distributionModel = new LogNormalDistributionModel(mean, stdev, clockParameterValues[3], lrcParametersInRealSpace);
	        
	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), //
	        		rateCategoryParameter, //
	                distributionModel, //
	                1, // 
	                false, // 
	                Double.NaN, //
	                true, //randomizeRates
	                false, // keepRates
					false // cacheRates
	                );

		} else if(this.clockModelIndex == 2) { // Exponential relaxed clock
		
			double numberOfBranches = 2 * (createTreeModel().getTaxonCount() - 1);
			Parameter rateCategoryParameter = new Parameter.Default(numberOfBranches);
			
			Parameter mean = new Parameter.Default(DistributionModelParser.MEAN, 1, clockParameterValues[4]);
	        ParametricDistributionModel distributionModel = new ExponentialDistributionModel(mean, clockParameterValues[5]);
			
//	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), rateCategoryParameter, 
//	                distributionModel, 1, false, Double.NaN);
			
	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), //
	        		rateCategoryParameter, //
	                distributionModel, //
	                1, // 
	                false, // 
	                Double.NaN, //
	                true, //randomizeRates
					false, // keepRates
					false // cacheRates
	                );
	        
		} else if(this.clockModelIndex == 3) { // Inverse Gaussian

			double numberOfBranches = 2 * (createTreeModel().getTaxonCount() - 1);
			Parameter rateCategoryParameter = new Parameter.Default(numberOfBranches);
			
			Parameter mean = new Parameter.Default(InverseGaussianDistributionModelParser.MEAN, 1, clockParameterValues[6]);
			Parameter stdev = new Parameter.Default(InverseGaussianDistributionModelParser.STDEV, 1, clockParameterValues[7]);
	        ParametricDistributionModel distributionModel = new InverseGaussianDistributionModel(
					mean, stdev, clockParameterValues[8], false);
     
	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), //
	        		rateCategoryParameter, //
	                distributionModel, //
	                1, // 
	                false, // 
	                Double.NaN, //
	                true, //randomizeRates
					false, // keepRates
					false // cacheRates
	                );
	        
		} else {

			System.out.println("Not yet implemented");

		}

		return branchRateModel;
	}// END: createBranchRateModel

	// ///////////////////////
	// ---SITE RATE MODEL---//
	// ///////////////////////

	public int siteRateModelIndex = 0;

	public String siteRateModelIdref = SiteModel.SITE_MODEL;

	public void resetSiteRateModelIdref() {
		this.siteRateModelIdref = SiteModel.SITE_MODEL;
	}
	
	public static String[] siteRateModels = { "No Model", //
			"Gamma Site Rate Model" //
	};

	public static String[] siteRateModelParameterNames = new String[] {
			"Gamma categories", // Gamma Site Rate Model
			"Alpha", // Gamma Site Rate Model
			"Invariant sites proportion" // Gamma Site Rate Model
	};

	public static int[][] siteRateModelParameterIndices = { {}, // NoModel
			{ 0, 1, 2 }, // GammaSiteRateModel
	};

	public double[] siteRateModelParameterValues = new double[] { 4.0, // GammaCategories
			0.5, // Alpha
			0.0 // Invariant sites proportion
	};

	public GammaSiteRateModel createSiteRateModel() {

		GammaSiteRateModel siteModel = null;
		String name = "siteModel";

		if (this.siteRateModelIndex == 0) { // no model

			siteModel = new GammaSiteRateModel(name);

		} else if (this.siteRateModelIndex == 1) { // GammaSiteRateModel

			siteModel = new GammaSiteRateModel(name,
					siteRateModelParameterValues[1],
					(int) siteRateModelParameterValues[0], siteRateModelParameterValues[2]);

		} else {

			System.out.println("Not yet implemented");

		}

		return siteModel;
	}// END: createGammaSiteRateModel

	// //////////////////////////
	// ---ANCESTRAL SEQUENCE---//
	// //////////////////////////
	
	public String ancestralSequenceString = null;

	public Sequence createAncestralSequence() {

		Sequence sequence = new Sequence(ancestralSequenceString);
		// sequence.appendSequenceString(ancestralSequenceString);

		return sequence;
	}

}// END: class

