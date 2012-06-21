package dr.app.bss;

import java.io.File;

import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;


public class BeagleSequenceSimulatorData {
	
	public static final String version = "1.0";
    public File treeFile = null;
    public TaxonList taxonList = new Taxa();
    public int replicateCount = 1000;
    public int popSizeCount = 0;
    public int substitutionModel = 0;
    
	/////////////////////////////
	//---SUBSTITUTION MODELS---//
	/////////////////////////////
    
	public static  String[] substitutionModels = { "HKY", //
			"Yang Codon Model" //
	};

	public static String[] substitutionParameterNames = new String[] { "Kappa value", // 
			"Omega value" //
	};

	public int[][] substitutionParameterIndices = { { 0 }, // HKY
			{ 0, 1 }, // Yang Codon Model
	};

	public double[] substitutionParameterValues = new double[] { 10.0, // Kappa value
			0.1, // Omega value
	};
	
	//////////////////////
	//---CLOCK MODELS---//
	//////////////////////
	
	// TODO: distributions for unrc

	public static  String[] clockModels = { "Strict Clock", //
	};
	
	public static String[] clockParameterNames = new String[] { "Clock rate", //
	};
	
	public int[][] clockParameterIndices = { { 0 }, // Strict Clock
	};

	public double[] clockParameterValues = new double[] { 1.2E-2, // clock rate
	};

	
	
	
	
	
	
	public BeagleSequenceSimulatorData() {
	}// END: Constructor
	
    public void setSubstitutionModel(int model) {
    	
    	substitutionModel = model;

    }//END: setSubstitutionModelModel
	
	
	
	
	
}//END: class

