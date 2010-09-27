package test.dr.integration;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.util.Utils;
import dr.evomodel.substmodel.NucModelType;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceException;
import test.dr.app.beauti.BeautiTesterConfig;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;



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
	
	private static final int decimal = 3;
	private static final String SEQUNECE_LEN = "1000";
	private static final String PRE_PATH = "C:\\Users\\dxie004\\workspace\\BEAST_MCMC\\examples\\Joseph\\test\\";
	private List<BEASTInputFile> inputFiles;
	
	
	public GTRParameterEstimationTest() {
		// Runtime.getRuntime().exec not work, have to run in Debug and run *.cmd manually each stage.
//		try{
//			splitTreeFiles();			
//		} catch (IOException ioe) {
//            System.err.println("Unable to read or write files: " + ioe.getMessage());
//        }
//		System.out.println(PRE_PATH + "seqgen_banch.cmd is generated"); 
//		
////	    try {
////	    	System.out.println(PRE_PATH + "seqgen_banch.cmd");
////	        Process p = Runtime.getRuntime().exec("cmd.exe " + PRE_PATH + "seqgen_banch.cmd");
////	        
////	        p.waitFor();
////	        System.out.println(p.exitValue());
////	    } catch (Exception err) {
////	    	System.err.println("Can not create seqgen_banch.cmd " + err.getMessage());
////		}      	    
////		System.out.println(PRE_PATH + "seqgen_banch.cmd is executed"); 
//		
//		createBeastInputFiles();
//		
//		System.out.println(PRE_PATH + "beast_banch.cmd is generated");
		
//	    try {
//	    	
//	        Process p = Runtime.getRuntime().exec("cmd.exe " + PRE_PATH + "beast_banch.cmd");
//	        
//	        p.waitFor();
//	        System.out.println(p.exitValue());
//	    } catch (Exception err) {
//	    	System.err.println("Can not create beast_banch.cmd " + err.getMessage());
//	    }      	    
//	    System.out.println(PRE_PATH + "beast_banch.cmd is executed"); 
		
		initTest();
		
		try{
			writeReport();
		} catch (IOException ioe) {
            System.err.println("Unable to write reprot file: " + ioe.getMessage());
        } catch (TraceException te) {
        	System.err.println("Problem to analyse log file: " + te.getMessage());
        }
		
		System.out.println("report is saved in " + PRE_PATH);
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
        		
        		b.setAc(roundDecimals(decimal, getRandomNum(0.05, 0.5)));
                b.setAg(roundDecimals(decimal, getRandomNum(0.3333, 3)));
                b.setAt(roundDecimals(decimal, getRandomNum(0.05, 0.5)));
                b.setCg(roundDecimals(decimal, getRandomNum(0.05, 0.5)));
                b.setGt(roundDecimals(decimal, getRandomNum(0.05, 0.5)));
                
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
        
//        btc.printlnScriptWriter("exit");
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
        btc.createScriptWriter(PRE_PATH + "beast_banch.cmd");
    	
        for (BEASTInputFile b : inputFiles) {
			beautiOptions = btc.createOptions();
	        btc.importFromFile(PRE_PATH + b.getFileNamePrefix() + ".nex", beautiOptions, false); 
	        
	        // assume only 1 partition model
	        PartitionSubstitutionModel model = beautiOptions.getPartitionSubstitutionModels().get(0);
	        
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
//	        beautiOptions.nodeHeightPrior = TreePrior.YULE;
//	        beautiOptions.clockType = ClockType.STRICT_CLOCK;
	        
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
	        beautiOptions.echoEvery = 10000;
	        beautiOptions.fileNameStem = PRE_PATH + b.getFileNamePrefix();
	        
	        
	        btc.setCOMMEND("java -jar beast.jar ");
	        btc.generate(PRE_PATH + b.getFileNamePrefix(), beautiOptions); // include btc.printlnScriptWriter( )       
	        
		}
//        btc.printlnScriptWriter("exit");
        btc.closeScriptWriter(); 
       
    }

    public void writeReport() throws IOException, TraceException {
    	FileWriter fileWriter = new FileWriter (PRE_PATH + "report.txt");
    	// columns 25
    	fileWriter.write("YULE+STRICT_CLOCK\tBIRTH_RATE_1\tBIRTH_RATE_2\tBIRTH_RATE_IN_RANGE\tBIRTH_RATE_DSS" + 
    			"\tAC_1\tAC_2\tAC_IN_RANGE\tAC_DSS" + "\tAG_1\tAG_2\tAG_IN_RANGE\tAG_DSS" + 
    			"\tAT_1\tAT_2\tAT_IN_RANGE\tAT_DSS" + "\tCG_1\tCG_2\tCG_IN_RANGE\tCG_DSS" + "\tGT_1\tGT_2\tGT_IN_RANGE\tGT_DSS\n"); 
		
    	String resultLine;
    	
    	for (BEASTInputFile b : inputFiles) {
    		resultLine = b.getFileNamePrefix();
    		
    		File logFile = new File(resultLine + ".log");
    		LogFileTraces traces = new LogFileTraces(resultLine + ".log", logFile);
            traces.loadTraces();
            traces.setBurnIn(100000);
//TODO: incorrect mean, stdv, ess, ...
            for (int i = 0; i < traces.getTraceCount(); i++) {
                traces.analyseTrace(i);
//            }
//    		
//            for (int i = 0; i < traces.getTraceCount(); i++) { 
            	
            	TraceCorrelation distribution = traces.getCorrelationStatistics(i);
//            	TraceDistribution distribution = traces.getDistributionStatistics(i);
            	
            	double hpdLower = distribution.getLowerHPD();
            	double hpdUpper = distribution.getUpperHPD();
            	
            	if (traces.getTraceName(i).equalsIgnoreCase("yule.birthRate")) {
            		
            		resultLine = resultLine + "\t" + Double.toString(BEASTInputFile.birthRate) + "\t" + Double.toString(distribution.getMean())
    									+ "\t" + Boolean.toString(isInRange(BEASTInputFile.birthRate, hpdLower, hpdUpper)) 
    									+ "\t" + Double.toString(distribution.getESS());
            		
            	} else if (traces.getTraceName(i).equalsIgnoreCase("ac")) {
            		
            		resultLine = resultLine + "\t" + Double.toString(b.getAc()) + "\t" + Double.toString(distribution.getMean())
					+ "\t" + Boolean.toString(isInRange(b.getAc(), hpdLower, hpdUpper)) 
					+ "\t" + Double.toString(distribution.getESS());
            		
            	} else if (traces.getTraceName(i).equalsIgnoreCase("ag")) {
            		
            		resultLine = resultLine + "\t" + Double.toString(b.getAg()) + "\t" + Double.toString(distribution.getMean())
					+ "\t" + Boolean.toString(isInRange(b.getAg(), hpdLower, hpdUpper)) 
					+ "\t" + Double.toString(distribution.getESS());
            		
            	} else if (traces.getTraceName(i).equalsIgnoreCase("at")) {
            		
            		resultLine = resultLine + "\t" + Double.toString(b.getAt()) + "\t" + Double.toString(distribution.getMean())
					+ "\t" + Boolean.toString(isInRange(b.getAt(), hpdLower, hpdUpper)) 
					+ "\t" + Double.toString(distribution.getESS());
            		
            	} else if (traces.getTraceName(i).equalsIgnoreCase("cg")) {
            		
            		resultLine = resultLine + "\t" + Double.toString(b.getCg()) + "\t" + Double.toString(distribution.getMean())
					+ "\t" + Boolean.toString(isInRange(b.getCg(), hpdLower, hpdUpper)) 
					+ "\t" + Double.toString(distribution.getESS());
            		
            	} else if (traces.getTraceName(i).equalsIgnoreCase("gt")) {
            		
            		resultLine = resultLine + "\t" + Double.toString(b.getGt()) + "\t" + Double.toString(distribution.getMean())
					+ "\t" + Boolean.toString(isInRange(b.getGt(), hpdLower, hpdUpper)) 
					+ "\t" + Double.toString(distribution.getESS());
            		
            	} 
            	
            	System.out.println("i = " + i + "   name = " + traces.getTraceName(i) + "  mean = " + distribution.getMean() + 
            			"   stdv = " + distribution.getStdErrorOfMean() + "   median = " + distribution.getMedian() +
            			"   min = " + distribution.getMinimum() + "   max = " + distribution.getMaximum() + 
            			"  ESS = " + distribution.getESS() + "   Lhpd = " + distribution.getLowerHPD() + "   Hhpd = " + distribution.getUpperHPD());
            }
              		
    		fileWriter.write(resultLine + "\n");
    	}
    	
    	fileWriter.close();
    }
    
    private boolean isInRange (double value, double min, double max) {
    	if (max > min) {
    		return (value > min && value < max);
    	} else {
    		return (value < min && value > max);
    	}
    }
    
    private void initTest() {
    	inputFiles = new ArrayList<BEASTInputFile>();
    	
    	BEASTInputFile b = new BEASTInputFile(); 
    			
		b.setAc(0.177);
        b.setAg(0.377);
        b.setAt(0.221);
        b.setCg(0.437);
        b.setGt(0.133);        
        b.setFileNamePrefix("STATE_0");
        
		inputFiles.add(b);
    	
    	
    }
    
    //Main method
    public static void main(String[] args) {

        new GTRParameterEstimationTest();
    }


 
}
