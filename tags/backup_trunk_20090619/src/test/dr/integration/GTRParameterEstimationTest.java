package test.dr.integration;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import test.dr.beauti.BeautiTesterConfig;
import dr.app.beauti.options.*;
import dr.app.util.Utils;


/**
 * GTR Parameter Estimation Tester.
 * Only available for Windows OS
 *
 * @author Walter Xie
 * @version 1.0
 * @since <pre>08/06/2009</pre>
 */
public class GTRParameterEstimationTest {

	
//    private static final String TREE_HEIGHT = CoalescentSimulator.ROOT_HEIGHT;
//    private static final String birthRateIndicator = "birthRateIndicator";
//    private static final String birthRate = "birthRate";
	
	private static final String SEQUNECE_LEN = "1000";
	private static final String PRE_PATH = "C:\\Users\\dxie004\\workspace\\BEAST_MCMC\\examples\\Joseph\\test\\";
	private List<BEASTInputFile> inputFiles;
	
	
	public GTRParameterEstimationTest() {
		
		try{
			splitTreeFiles();			
		} catch (IOException ioe) {
            System.err.println("Unable to read or write files: " + ioe.getMessage());
        }
		System.out.println(PRE_PATH + "seqgen_banch.cmd is generated"); 
		
//	    try {
//	        Process p = Runtime.getRuntime().exec(PRE_PATH + "seqgen_banch.cmd");
//	        
//	        p.waitFor();
//	        System.out.println(p.exitValue());
//	    } catch (Exception err) {
//	    	System.err.println("Can not create seqgen_banch.cmd " + err.getMessage());
//		}      	    
		System.out.println(PRE_PATH + "seqgen_banch.cmd is executed"); 
		
//		createBeastXMLFiles();
	}
    
	
	public void splitTreeFiles() throws IOException {
		inputFiles = new ArrayList<BEASTInputFile>();
		
		FileReader fileReader = new FileReader(PRE_PATH + "n.trees");
		BufferedReader reader = new BufferedReader(fileReader);
		
		FileWriter fileWriter; 
//		FileWriter cmdWriter = new FileWriter (PRE_PATH + "seqgen_banch.cmd");
		
		BeautiTesterConfig btc = new BeautiTesterConfig();
        btc.createScriptWriter(PRE_PATH + "seqgen_banch.cmd");
                       
        String line = Utils.nextNonCommentLine(reader);
        
        while (line != null) {
            
        	if (line.startsWith("tree")) { 
        		BEASTInputFile b = new BEASTInputFile(); 
        		// read 1 tree each line 
        		String[] values = line.split(" ");
        		// write this tree in a file
        		fileWriter = new FileWriter (PRE_PATH + values[1] + ".tree"); // STATE_1000
        		fileWriter.write(values[values.length - 1]); // tree
        		fileWriter.close();
        		
        		b.setAc(roundDecimals(5, getRandomNum(0.05, 0.5)));
                b.setAg(roundDecimals(5, getRandomNum(0.3333, 3)));
                b.setAt(roundDecimals(5, getRandomNum(0.05, 0.5)));
                b.setCg(roundDecimals(5, getRandomNum(0.05, 0.5)));
                b.setGt(roundDecimals(5, getRandomNum(0.05, 0.5)));
                
                b.setFileNamePrefix(values[1].trim());
                
        		inputFiles.add(b);
        		
        		String comLine = "C:\\Users\\dxie004\\Documents\\Seq-Gen.v1.3.2\\seq-gen -mGTR -fe -on -l" + SEQUNECE_LEN + " -r" + 
        				Double.toString(b.getAc()) + "," + Double.toString(b.getAg()) + "," + Double.toString(b.getAt()) + "," + 
        				Double.toString(b.getCg()) + "," + Double.toString(BEASTInputFile.ct) + "," + Double.toString(b.getGt()) + 
        				" < " + b.getFileNamePrefix() + ".tree > " + b.getFileNamePrefix() + ".nex";
// C:\Users\dxie004\Documents\Seq-Gen.v1.3.2\seq-gen -mGTR -fe -on -l1000 -r0.1,2.1,0.2,0.08,1,0.36 < STATE_0.tree > STATE_0.nex
          		
        		btc.printlnScriptWriter(comLine);    
//        		cmdWriter.write(comLine);
//        		cmdWriter.write("\n\n");
//        		System.out.println(comLine); 
        	}

            line = Utils.nextNonCommentLine(reader);
        }        
        
        btc.closeScriptWriter();
//        cmdWriter.close();
        fileReader.close();
	}
	
	private double getRandomNum(double min, double max) { // range
		Random r = new Random(); 
		if (max > min) {
			return (max-min)*r.nextDouble() + min;
		} else {
			return (min-max)*r.nextDouble() + max;
		}
	}
	
	private double roundDecimals(int decim, double d) {
		String tmp = "#.";
		if (decim < 1) {
			tmp = "#.#"; 
		}
		
		for (int i = 0; i < decim; i++) {
			tmp = tmp + "#";
		}
		
    	DecimalFormat twoDForm = new DecimalFormat(tmp);
    	return Double.valueOf(twoDForm.format(d));
	}

    public void createBeastInputFiles() {
    	BeautiOptions beautiOptions;
        BeautiTesterConfig btc = new BeautiTesterConfig();
        btc.createScriptWriter(PRE_PATH + "beast_banch.bat");
    	
        for (BEASTInputFile b : inputFiles) {
			beautiOptions = btc.createOptions();
	        btc.importFromFile(PRE_PATH + b.getFileNamePrefix() + ".nex", beautiOptions, false); 
	        
	        // assume only 1 partition model
	        PartitionModel model = beautiOptions.getPartitionModels().get(0);
	        
	        // set substitution model 
	        model.setNucSubstitutionModel(NucModelType.GTR);
//	        model.setNucSubstitutionModel(NucModelType.HKY);
	        model.setCodonHeteroPattern(null);
	        model.setUnlinkedSubstitutionModel(false);
	        model.setUnlinkedHeterogeneityModel(false);
	        
	        model.setGammaHetero(false);
	        model.setGammaCategories(4);
	        model.setInvarHetero(false);
	        
	        // set tree prior
	        beautiOptions.nodeHeightPrior = TreePrior.YULE;
	        beautiOptions.clockType = ClockType.STRICT_CLOCK;
	        
	        // change value of parameters
//	        Parameter para = beautiOptions.getParameter("yule.birthRate");
//	        para.initial = 50;
//	        // for multi-partition models, use beautiOptions.getParameter(name, model);
//	        
//	        // remove operators
//	        Operator op = beautiOptions.getOperator("yule.birthRate");
//	        op.inUse = false;	        
//	        op = model.getOperator("kappa");
//	        op.inUse = false;
//	        op = model.getOperator("frequencies");
//	        op.inUse = false;
	        
	        // set MCMC
	        beautiOptions.chainLength = 1000000;
	        beautiOptions.logEvery = 100;
	        beautiOptions.fileNameStem = PRE_PATH + b.getFileNamePrefix();
	        
	        
	        btc.setCOMMEND("java -jar beast.jar ");
	        btc.generate(PRE_PATH + b.getFileNamePrefix(), beautiOptions); // include btc.printlnScriptWriter( )       
	        
		}
   
        btc.closeScriptWriter(); 
       
    }

    public void writeReport() {
    	
    }
    
    
    
    //Main method
    public static void main(String[] args) {

        new GTRParameterEstimationTest();
    }


 
}
