package dr.app.bss;

import java.io.File;
import java.io.Serializable;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.GTR;
import dr.app.beagle.evomodel.substmodel.GY94CodonModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.TN93;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

@SuppressWarnings("serial")
public class PartitionData implements Serializable {

	public PartitionData() {
	}// END: Constructor
	
	public int from = 1;
	public int to = 1000;
	public int every = 1;
	
	// //////////////////
	// ---TREE MODEL---//
	// //////////////////
	
	public File treeFile = null;
	public transient String treeModelIdref =  TreeModel.TREE_MODEL;

	public TreeModel createTreeModel() {
		TreeModel treeModel = Utils.importTreeFromFile(treeFile);
		return treeModel;
	}
	
	// /////////////////
	// ---DATA TYPE---//
	// /////////////////
	
	public int dataTypeIndex = 0;

	public static String[] dataTypes = { "Nucleotide", //
			"Codon" //
	};
	
	public DataType createDataType() {

		DataType dataType = null;

		if (this.dataTypeIndex == 0) { // Nucleotide

			dataType = Nucleotides.INSTANCE;

		} else if (this.dataTypeIndex == 1) { // Codon

			dataType = Codons.UNIVERSAL;

		} else {

			System.out.println("Not yet implemented");

		}

		return dataType;
	}// END: createDataType
	
	// ///////////////////////////
	// ---SUBSTITUTION MODELS---//
	// ///////////////////////////
	
	public int substitutionModelIndex = 0;

	public transient String substitutionModelIdref = createSubstitutionModelIdref();
	
	private String createSubstitutionModelIdref() {

		String substitutionModelIdref = "";

		switch (this.substitutionModelIndex) {

		case 0: // HKY

			substitutionModelIdref = PartitionData.substitutionModels[0]
					.toLowerCase();
			break;

		case 1: // GTR

			substitutionModelIdref = PartitionData.substitutionModels[1]
					.toLowerCase();
			break;

		case 2: // TNF93
			substitutionModelIdref = PartitionData.substitutionModels[2].toLowerCase();
			break;

		case 3: // Yang Codon Model

			substitutionModelIdref = PartitionData.substitutionModels[3].replaceAll(" +", ".").toLowerCase();
			break;

		}// END: switch

		return substitutionModelIdref;
	}

	public static String[] substitutionModels = { "HKY", //
		    "GTR", //
		    "TN93", //
			"Yang Codon Model" //
	};

	public static String[] substitutionParameterNames = new String[] {
			"Kappa value", // HKY
			"AC", // GTR
			"AG", // GTR
			"AT", // GTR
			"CG", // GTR
			"CT", // GTR
			"GT", // GTR
			"Kappa 1", // TN93
			"Kappa 2", // TN93
			"Omega value", // Yang Codon Model
			"Kappa value" // Yang Codon Model
	};

	public int[][] substitutionParameterIndices = { { 0 }, // HKY
			{ 1, 2, 3, 4, 5, 6 }, // GTR
			{ 7, 8 }, // TN93
			{ 9, 10 }, // Yang Codon Model
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
			1.0 // Kappa value
	};

	public BranchModel createBranchModel() {

		BranchModel branchModel = null;

		if (this.substitutionModelIndex == 0) { // HKY

			Parameter kappa = new Parameter.Default(1, substitutionParameterValues[0]);
			
			FrequencyModel frequencyModel = this.createFrequencyModel();

			HKY hky = new HKY(kappa, frequencyModel);
			
			branchModel = new HomogeneousBranchModel(hky);

		} else if (this.substitutionModelIndex == 1) { // GTR

			Parameter ac = new Parameter.Default(1, substitutionParameterValues[1]);
			Parameter ag = new Parameter.Default(1, substitutionParameterValues[2]);
			Parameter at = new Parameter.Default(1, substitutionParameterValues[3]);
			Parameter cg = new Parameter.Default(1, substitutionParameterValues[4]);
			Parameter ct = new Parameter.Default(1, substitutionParameterValues[5]);
			Parameter gt = new Parameter.Default(1, substitutionParameterValues[6]);
			
			FrequencyModel frequencyModel = this.createFrequencyModel();
			
			GTR gtr = new GTR(ac, ag, at, cg, ct, gt, frequencyModel);

			branchModel = new HomogeneousBranchModel(gtr);
			
		} else if (this.substitutionModelIndex == 2) { // TN93

			Parameter kappa1 = new Parameter.Default(1, substitutionParameterValues[7]);
			Parameter kappa2 = new Parameter.Default(1, substitutionParameterValues[8]);
			
			FrequencyModel frequencyModel = this.createFrequencyModel();

			TN93 tn93 = new TN93(kappa1, kappa2, frequencyModel);
			
			branchModel = new HomogeneousBranchModel(tn93);
			
		} else if (this.substitutionModelIndex == 3) { // Yang Codon Model

            FrequencyModel frequencyModel = this.createFrequencyModel();

			Parameter kappa = new Parameter.Default(1, substitutionParameterValues[9]);
			Parameter omega = new Parameter.Default(1, substitutionParameterValues[10]);
			
			GY94CodonModel yangCodonModel = new GY94CodonModel(Codons.UNIVERSAL, omega, kappa, frequencyModel);
			
			branchModel = new HomogeneousBranchModel(yangCodonModel);

		} else if (this.substitutionModelIndex == 4) { 
			
			System.out.println("Not yet implemented");
			
		}

		return branchModel;
	}// END: createBranchSubstitutionModel

	// ////////////////////
	// ---CLOCK MODELS---//
	// ////////////////////

	public int clockModelIndex = 0;

	public transient String clockModelIdref = BranchRateModel.BRANCH_RATES;
	
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

		if (this.clockModelIndex == 0) { // Strict Clock

			Parameter rateParameter = new Parameter.Default(1, clockParameterValues[0]);
			branchRateModel = new StrictClockBranchRates(rateParameter);

		} else if (this.clockModelIndex == 1) {

			System.out.println("Not yet implemented");

		}

		return branchRateModel;
	}// END: createBranchRateModel
	
	// ////////////////////////
	// ---FREQUENCY MODELS---//
	// ////////////////////////

	public transient String frequencyModelIdref = "freqModel";
	
	public int frequencyModelIndex = 0;

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

		if (this.frequencyModelIndex == 0) { // Nucleotidefrequencies

			Parameter freqs = new Parameter.Default(new double[] {
					frequencyParameterValues[0], frequencyParameterValues[1],
					frequencyParameterValues[2], frequencyParameterValues[3] });

			frequencyModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);

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

			frequencyModel = new FrequencyModel(Codons.UNIVERSAL, freqs);

		} else if (this.frequencyModelIndex == 2) {
			
			System.out.println("Not yet implemented");
			
		}

		return frequencyModel;
	}// END: createFrequencyModel

	// ////////////////////////
	// ---SITE RATE MODELS---//
	// ////////////////////////

	public int siteRateModelIndex = 0;

	public transient String siteRateModelIdref = SiteModel.SITE_MODEL;
	
	public static String[] siteRateModels = { "No model", //
			"Gamma Site Rate Model", //
	};

	public static String[] siteRateModelParameterNames = new String[] {
			"Gamma categories", //
			"alpha", //
	};

	public int[][] siteRateModelParameterIndices = { {}, // nomodel
			{ 0, 1 }, // GammaSiteRateModel
	};

	public double[] siteRateModelParameterValues = new double[] { 4.0, // GammaCategories
			0.5, // alpha
	};

	public GammaSiteRateModel createSiteRateModel() {

		GammaSiteRateModel siteModel = null;
		String name = "siteModel";

		if (this.siteRateModelIndex == 0) { // no model

			siteModel = new GammaSiteRateModel("siteModel");

		} else if (this.siteRateModelIndex == 1) { // GammaSiteRateModel

			siteModel = new GammaSiteRateModel(name, siteRateModelParameterValues[1],
					(int) siteRateModelParameterValues[0]);

		} else if (this.siteRateModelIndex == 2) {

			System.out.println("Not yet implemented");

		}

		return siteModel;
	}// END: createGammaSiteRateModel
	
}// END: class

