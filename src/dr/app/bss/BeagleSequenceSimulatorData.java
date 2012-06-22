package dr.app.bss;

import java.io.File;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.HomogenousBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

public class BeagleSequenceSimulatorData {

	public static final String version = "1.0";
	public File treeFile = null;
	public TaxonList taxonList = new Taxa();
	public int replicateCount = 1000;
	public TreeModel treeModel;
	
	// ///////////////////////////
	// ---SUBSTITUTION MODELS---//
	// ///////////////////////////

	public int substitutionModel = 0;

	public static String[] substitutionModels = { "HKY", //
			"Yang Codon Model" //
	};

	public static String[] substitutionParameterNames = new String[] {
			"Kappa value", //
			"Omega value" //
	};

	public int[][] substitutionParameterIndices = { { 0 }, // HKY
			{ 0, 1 }, // Yang Codon Model
	};

	public double[] substitutionParameterValues = new double[] { 10.0, // Kappa-value
			0.1, // Omega value
	};

	public void setSubstitutionModel(int model) {
		substitutionModel = model;
	}// END: setSubstitutionModelModel
	
	public BranchSubstitutionModel createBranchSubstitutionModel() {

		BranchSubstitutionModel substitutionModel = null;

		if (this.substitutionModel == 0) {// HKY

			FrequencyModel frequencyModel = this.createFrequencyModel();
			Parameter kappa = new Parameter.Default(1, substitutionParameterValues[0]);
			HKY hky = new HKY(kappa, frequencyModel);
			substitutionModel = new HomogenousBranchSubstitutionModel(hky,
					frequencyModel);

		} else if (this.substitutionModel == 0) {

			System.out.println("Not yet implemented");

		}

		return substitutionModel;
	}// END: createBranchSubstitutionModel

	// ////////////////////
	// ---CLOCK MODELS---//
	// ////////////////////

	public int clockModel = 0;

	public static String[] clockModels = { "Strict Clock", //
	};

	public static String[] clockParameterNames = new String[] { "Clock rate", //
	};

	public int[][] clockParameterIndices = { { 0 }, // StrictClock
	};

	public double[] clockParameterValues = new double[] { 1.2E-2, // clockrate
	};

	public BranchRateModel createBranchRateModel() {

		BranchRateModel branchRateModel = null;

		if (this.clockModel == 0) {// Strict Clock

			Parameter rateParameter = new Parameter.Default(1, clockParameterValues[0]);
			branchRateModel = new StrictClockBranchRates(rateParameter);

		} else if (this.clockModel == 1) {

			System.out.println("Not yet implemented");

		}

		return branchRateModel;
	}// END: createBranchRateModel
	
	// ////////////////////////
	// ---FREQUENCY MODELS---//
	// ////////////////////////

	public int frequencyModel = 0;

	public static String[] frequencyModels = { "Nucleotide frequencies", //
		"Codon frequencies"
	};

	public static String[] frequencyParameterNames = new String[] {
			"Nucleotide frequencies 1", //
			"Nucleotide frequencies 2", //
			"Nucleotide frequencies 3", //
			"Nucleotide frequencies 4", //
			"Codon frequencies 1", //
			"Codon frequencies 2", //
			"Codon frequencies 3", //
			"Codon frequencies 4", //
			"Codon frequencies 5", //
			"Codon frequencies 6", //
			"Codon frequencies 7", //
			"Codon frequencies 8", //
			"Codon frequencies 9", //
			"Codon frequencies 10", //
			"Codon frequencies 11", //
			"Codon frequencies 12", //
			"Codon frequencies 13", //
			"Codon frequencies 14", //
			"Codon frequencies 15", //
			"Codon frequencies 16", //
			"Codon frequencies 17", //
			"Codon frequencies 18", //
			"Codon frequencies 19", //
			"Codon frequencies 20", //
			"Codon frequencies 21", //
			"Codon frequencies 22", //
			"Codon frequencies 23", //
			"Codon frequencies 24", //
			"Codon frequencies 25", //
			"Codon frequencies 26", //
			"Codon frequencies 27", //
			"Codon frequencies 28", //
			"Codon frequencies 29", //
			"Codon frequencies 30", //
			"Codon frequencies 31", //
			"Codon frequencies 32", //
			"Codon frequencies 33", //
			"Codon frequencies 34", //
			"Codon frequencies 35", //
			"Codon frequencies 36", //
			"Codon frequencies 37", //
			"Codon frequencies 38", //
			"Codon frequencies 39", //
			"Codon frequencies 40", //
			"Codon frequencies 41", //
			"Codon frequencies 42", //
			"Codon frequencies 43", //
			"Codon frequencies 44", //
			"Codon frequencies 45", //
			"Codon frequencies 46", //
			"Codon frequencies 47", //
			"Codon frequencies 48", //
			"Codon frequencies 49", //
			"Codon frequencies 50", //
			"Codon frequencies 51", //
			"Codon frequencies 52", //
			"Codon frequencies 53", //
			"Codon frequencies 54", //
			"Codon frequencies 55", //
			"Codon frequencies 56", //
			"Codon frequencies 57", //
			"Codon frequencies 58", //
			"Codon frequencies 59", //
			"Codon frequencies 60", //
			"Codon frequencies 61", //
	};

	public int[][] frequencyParameterIndices = { { 0, 1, 2, 3 }, // Nucleotidefrequencies
			{ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, //
					14, 15, 16, 17, 18, 19, 20, 21, 22, 23, //
					24, 25, 26, 27, 28, 29, 30, 31, 32, 33, //
					34, 35, 36, 37, 38, 39, 40, 41, 42, 43, //
					44, 45, 46, 47, 48, 49, 50, 51, 52, 53, //
					54, 55, 56, 57, 58, 59, 60, 61, 62, 63, //
					64 } // Codonfrequencies
	};

	public double[] frequencyParameterValues = new double[] { 0.25, // nucleotidefrequencies1
			0.25, // nucleotidefrequencies2
			0.25, // nucleotidefrequencies3
			0.25, // nucleotidefrequencies4
			0.0163936, // codonfrequencies1
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
			0.01639344 // codonfrequencies61
	};

	public FrequencyModel createFrequencyModel() {

		FrequencyModel frequencyModel = null;

		if (this.frequencyModel == 0) { // Nucleotidefrequencies

			Parameter freqs = new Parameter.Default(new double[] {
					frequencyParameterValues[0], frequencyParameterValues[1],
					frequencyParameterValues[2], frequencyParameterValues[3] });
			frequencyModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);

		} else if (this.frequencyModel == 1) {

			System.out.println("Not yet implemented");

		}

		return frequencyModel;
	}// END: createFrequencyModel
	
	// ////////////////////////
	// ---SITE RATE MODELS---//
	// ////////////////////////

	public int siteModel = 0;

	public static String[] siteModels = { "No model", //
			"Gamma Site Rate Model", //
	};

	public static String[] siteParameterNames = new String[] {
			"Gamma categories", //
			"alpha", //
	};

	public int[][] siteParameterIndices = { {}, // nomodel
			{ 0, 1 }, // GammaSiteRateModel
	};

	public double[] siteParameterValues = new double[] { 4.0, // GammaCategories
			0.5, // alpha
	};

	public GammaSiteRateModel createSiteRateModel() {

		GammaSiteRateModel siteModel = null;
		String name = "siteModel";

		if (this.siteModel == 0) { // no model

			siteModel = new GammaSiteRateModel("siteModel");

		} else if (this.siteModel == 1) { // GammaSiteRateModel

			siteModel = new GammaSiteRateModel(name, siteParameterValues[1],
					(int) siteParameterValues[0]);

		} else if (this.siteModel == 2) {

			System.out.println("Not yet implemented");

		}

		return siteModel;
	}// END: createGammaSiteRateModel
	
	public BeagleSequenceSimulatorData() {
	}// END: Constructor

}// END: class

