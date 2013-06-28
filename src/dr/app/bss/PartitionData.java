package dr.app.bss;

import java.io.Serializable;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.GTR;
import dr.app.beagle.evomodel.substmodel.GY94CodonModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.TN93;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.SiteModel;
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
	public int to = 1000;
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

	// /////////////////////////
	// ---DEMOGRAPHIC MODEL---//
	// /////////////////////////

	public int demographicModelIndex = 0;
	
	public String demographicModelIdref = Utils.DEMOGRAPHIC_MODEL;

	public void resetDemographicModelIdref() {
		this.demographicModelIdref = Utils.DEMOGRAPHIC_MODEL;
	}
	
    public static String[] demographicModels = {
    	"No Model (user-specified tree)",
    	"Constant Population",
        "Exponential Growth (Growth Rate)",
        "Exponential Growth (Doubling Time)"
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
			"Doubling Time" // Exponential Growth (Doubling Time)
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
			{ 3, 4 }// Exponential Growth (Doubling Time)
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
			10.0 // Doubling Time
			
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
	
//	public Taxa taxa = null;
	
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
			
		} else if( (this.demographicModelIndex > 0 && this.demographicModelIndex <= 3) && this.record.isTreeSet()) {
			
			Taxa taxa = new Taxa(this.record.getTree().asList()); 
			CoalescentSimulator topologySimulator = new CoalescentSimulator();
			treeModel = new TreeModel(topologySimulator.simulateTree(taxa, createDemographicFunction()));			
			
		} else if((this.demographicModelIndex > 0 && this.demographicModelIndex <= 3) && this.record.isTaxaSet()) {
			
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

//	public int dataTypeIndex = 0;
//
//	public static String[] dataTypes = { "Nucleotide", //
//			"Codon" //
//	};
//
//	public DataType createDataType() {
//
//		DataType dataType = null;
//
//		if (this.dataTypeIndex == 0) { // Nucleotide
//
//			dataType = Nucleotides.INSTANCE;
//
//		} else if (this.dataTypeIndex == 1) { // Codon
//
//			dataType = Codons.UNIVERSAL;
//
//		} else {
//
//			System.out.println("Not yet implemented");
//
//		}
//
//		return dataType;
//	}// END: createDataType

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
			"GY94CodonModel" //
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

	public static int[][] substitutionParameterIndices = { { 0 }, // HKY
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

		} else {

			System.out.println("Not yet implemented");

		}

		return branchModel;
	}// END: createBranchSubstitutionModel

	// ////////////////////
	// ---CLOCK MODELS---//
	// ////////////////////
	
	public int clockModelIndex = 0;

	public String clockModelIdref = BranchRateModel.BRANCH_RATES;

	public void resetClockModelIdref() {
		this.clockModelIdref = BranchRateModel.BRANCH_RATES;
	}
	
	public static String[] clockModels = { "Strict Clock", //
			"Lognormal relaxed clock (Uncorrelated)", //
			"Exponential relaxed clock (Uncorrelated)", //
			"Inverse Gaussian relaxed clock" //
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

	public double[] clockParameterValues = new double[] { 1.2E-2, // clockrate
			1.0, // ucld.mean
			2.0, // ucld.stdev
			0.0, // ucld.offset
			1.0, // uced.mean
            0.0, // uced.offset
			0.0, // ig.mean
			1.0, // ig.stdev
			0.0 // ig.offset
	};

	public BranchRateModel createClockRateModel() {

		BranchRateModel branchRateModel = null;

		if (this.clockModelIndex == 0) { // Strict Clock

			Parameter rateParameter = new Parameter.Default(1, clockParameterValues[0]);
			
			branchRateModel = new StrictClockBranchRates(rateParameter);

		} else if (this.clockModelIndex == 1) {// Lognormal relaxed clock

			double numberOfBranches = 2 * (createTreeModel().getTaxonCount() - 1);
			Parameter rateCategoryParameter = new Parameter.Default(numberOfBranches);
			
			Parameter mean = new Parameter.Default(LogNormalDistributionModelParser.MEAN, 1, clockParameterValues[1]);
			Parameter stdev = new Parameter.Default(LogNormalDistributionModelParser.STDEV, 1, clockParameterValues[2]);
	        ParametricDistributionModel distributionModel = new LogNormalDistributionModel(mean, stdev, clockParameterValues[3], true, true);
	        
	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), rateCategoryParameter, 
	                distributionModel, 1, false, Double.NaN);

		} else if(this.clockModelIndex == 2) { // Exponential relaxed clock
		
			double numberOfBranches = 2 * (createTreeModel().getTaxonCount() - 1);
			Parameter rateCategoryParameter = new Parameter.Default(numberOfBranches);
			
			Parameter mean = new Parameter.Default(DistributionModelParser.MEAN, 1, clockParameterValues[4]);
	        ParametricDistributionModel distributionModel = new ExponentialDistributionModel(mean, clockParameterValues[5]);
			
	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), rateCategoryParameter, 
	                distributionModel, 1, false, Double.NaN);
			
		} else if(this.clockModelIndex == 3) { // Inverse Gaussian

			double numberOfBranches = 2 * (createTreeModel().getTaxonCount() - 1);
			Parameter rateCategoryParameter = new Parameter.Default(numberOfBranches);
			
			Parameter mean = new Parameter.Default(InverseGaussianDistributionModelParser.MEAN, 1, clockParameterValues[6]);
			Parameter stdev = new Parameter.Default(InverseGaussianDistributionModelParser.STDEV, 1, clockParameterValues[7]);
	        ParametricDistributionModel distributionModel = new InverseGaussianDistributionModel(
					mean, stdev, clockParameterValues[8], false);
     
	        branchRateModel = new DiscretizedBranchRates(createTreeModel(), rateCategoryParameter, 
	                distributionModel, 1, false, Double.NaN);
	        
		} else {

			System.out.println("Not yet implemented");

		}

		return branchRateModel;
	}// END: createBranchRateModel

	// ////////////////////////
	// ---FREQUENCY MODELS---//
	// ////////////////////////

	public String frequencyModelIdref = Utils.FREQUENCY_MODEL;

	public void resetFrequencyModelIdref() {
		this.frequencyModelIdref = Utils.FREQUENCY_MODEL;
	}
	
	public int frequencyModelIndex = 0;

	public static String[] frequencyModels = { "Nucleotide frequencies", //
			"Codon frequencies" };

	public static String[] frequencyParameterNames = new String[] {
			"A frequency", // Nucleotide frequencies
			"C frequency", // Nucleotide frequencies
			"G frequency", // Nucleotide frequencies
			"T frequency", // Nucleotide frequencies
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

	public double[] frequencyParameterValues = new double[] { 0.25, // A frequency
			0.25, // C frequency
			0.25, // G frequency
			0.25, // T frequency
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

		} else {

			System.out.println("Not yet implemented");

		}

		return frequencyModel;
	}// END: createFrequencyModel

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
//		sequence.appendSequenceString(ancestralSequenceString);
		
		return sequence;
	}
	
	
}// END: class

