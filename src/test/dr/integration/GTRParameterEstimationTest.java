package test.dr.integration;

import java.io.*;
import java.util.Arrays;

import jebl.math.Random;



import test.dr.beauti.BeautiTesterConfig;
import dr.app.beauti.options.*;
import dr.app.util.Arguments;
import dr.app.util.Utils;


/**
 * GTR Parameter Estimation Tester.
 *
 * @author Walter Xie
 * @version 1.0
 * @since <pre>08/06/2009</pre>
 */
public class GTRParameterEstimationTest {

	
//    private static final String TREE_HEIGHT = CoalescentSimulator.ROOT_HEIGHT;
//    private static final String birthRateIndicator = "birthRateIndicator";
//    private static final String birthRate = "birthRate";

	private static final String PRE_PATH = "examples/Joseph/test/";
	
	public GTRParameterEstimationTest() {
		
		try{
			splitTreeFiles();			
		} catch (IOException ioe) {
            System.err.println("Unable to read or write files: " + ioe.getMessage());
        }
		
//		createBeastXMLFiles();
	}
    
	
	public void splitTreeFiles() throws IOException {
		
		FileReader fileReader = new FileReader(PRE_PATH + "n.trees");
		BufferedReader reader = new BufferedReader(fileReader);
		
		FileWriter fileWriter; 
		
		BeautiTesterConfig btc = new BeautiTesterConfig();
        btc.createScriptWriter(PRE_PATH + "seqgen_banch.bat");
        
        double ac;
        double ag;
        double at;
        double cg;
        double ct = 1; // fixed
        double gt;
        
        String line = Utils.nextNonCommentLine(reader);
        
        while (line != null) {
            
        	if (line.startsWith("tree")) { 
        		// read 1 tree each line 
        		String[] values = line.split(" ");
        		// write this tree in a file
        		fileWriter = new FileWriter (PRE_PATH + values[1] + ".tree"); // STATE_1000
        		fileWriter.write(values[values.length - 1]); // tree
        		fileWriter.close();
        		
        		ac = getRandomNum(0.05, 0.5);
                ag = getRandomNum(0.3333, 3);
                at = getRandomNum(0.05, 0.5);
                cg = getRandomNum(0.05, 0.5);
                gt = getRandomNum(0.05, 0.5);
        		
        		btc.printlnScriptWriter("C:\\Users\\dxie004\\Documents\\Seq-Gen.v1.3.2\\seq-Gen.v1.3.2\\seq-gen -mGTR -fe -on -l1000 -r" + 
        				Double.toString(ac) + "," + Double.toString(ag) + "," + Double.toString(at) + "," + 
        				Double.toString(cg) + "," + Double.toString(ct) + "," + Double.toString(gt) + 
        				" < " + values[1] + ".tree > " + values[1] + ".nex");
// C:\Users\dxie004\Documents\Seq-Gen.v1.3.2\seq-Gen.v1.3.2\seq-gen -mGTR -fe -on -l1000 -r0.1,2.1,0.2,0.08,1,0.36 < STATE_0.tree > STATE_0.nex
        	}

            line = Utils.nextNonCommentLine(reader);
        }        
        
        btc.closeScriptWriter();
        fileReader.close();
	}
	
	private double getRandomNum(double min, double max) { // range
		if (max > min) {
			return (max-min)*Random.nextDouble() + min;
		} else {
			return (min-max)*Random.nextDouble() + max;
		}
	}

    public void createBeastXMLFiles() {
    	BeautiOptions beautiOptions;
        BeautiTesterConfig btc = new BeautiTesterConfig();
        btc.createScriptWriter(PRE_PATH + "beast_banch.bat");
    	
        for (int i = 0; i < 3; i++) {
			beautiOptions = btc.createOptions();
	        btc.importFromFile(PRE_PATH + "n.nex", beautiOptions, false); 
	        
	        // assume only 1 partition model
	        PartitionModel model = beautiOptions.getPartitionModels().get(0);
	        
	        // set substitution model 
//	        model.setNucSubstitutionModel(NucModelType.GTR);
	        model.setNucSubstitutionModel(NucModelType.HKY);
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
	        Parameter para = beautiOptions.getParameter("yule.birthRate");
	        para.initial = 50;
	        // for multi-partition models, use beautiOptions.getParameter(name, model);
	        
	        // remove operators
	        Operator op = beautiOptions.getOperator("yule.birthRate");
	        op.inUse = false;	        
	        op = model.getOperator("kappa");
	        op.inUse = false;
	        op = model.getOperator("frequencies");
	        op.inUse = false;
	        
	        // set MCMC
	        beautiOptions.chainLength = 2000000;;
	        
	        btc.setCOMMEND("java -jar beast.jar ");
	        btc.generate(PRE_PATH + "gtr_" + Integer.toString(i), beautiOptions);        
	        
		}
   
        btc.closeScriptWriter(); 
       
    }

    //Main method
    public static void main(String[] args) {

        new GTRParameterEstimationTest();
    }

//    private void randomLocalYuleTester(TreeModel treeModel, Parameter I, Parameter b, OperatorSchedule schedule) {
//
//        MCMC mcmc = new MCMC("mcmc1");
//        MCMCOptions options = new MCMCOptions();
//        options.setChainLength(1000000);
//        options.setUseCoercion(true);
//        options.setPreBurnin(100);
//        options.setTemperature(1.0);
//        options.setFullEvaluationCount(2000);
//
//        TreelengthStatistic tls = new TreelengthStatistic(TL, treeModel);
//        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);
//
//        Parameter m = new Parameter.Default("m", 1.0, 0.0, Double.MAX_VALUE);
//
//        SpeciationModel speciationModel = new RandomLocalYuleModel(b, I, m, false, Units.Type.YEARS, 4);
//
//        Likelihood likelihood = new SpeciationLikelihood(treeModel, speciationModel, "randomYule.like");
//
//        ArrayLogFormatter formatter = new ArrayLogFormatter(false);
//
//        MCLogger[] loggers = new MCLogger[2];
//        loggers[0] = new MCLogger(formatter, 100, false);
//        loggers[0].add(likelihood);
//        loggers[0].add(rootHeight);
//        loggers[0].add(tls);
//        loggers[0].add(I);
//
//        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 100000, false);
//        loggers[1].add(likelihood);
//        loggers[1].add(rootHeight);
//        loggers[1].add(tls);
//        loggers[1].add(I);
//
//        mcmc.setShowOperatorAnalysis(true);
//
//        mcmc.init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers);
//
//        mcmc.run();
//
//        List<Trace> traces = formatter.getTraces();
//        ArrayTraceList traceList = new ArrayTraceList("yuleModelTest", traces, 0);
//
//        for (int i = 1; i < traces.size(); i++) {
//            traceList.analyseTrace(i);
//        }
//
//        TraceCorrelation tlStats =
//                traceList.getCorrelationStatistics(traceList.getTraceIndex("root." + birthRateIndicator));
//
//        System.out.println("mean = " + tlStats.getMean());
//        System.out.println("expected mean = 0.5");
//
//        assertExpectation("root." + birthRateIndicator, tlStats, 0.5);
//    }

 
}
