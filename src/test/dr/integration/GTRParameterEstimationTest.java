package test.dr.integration;

import test.dr.beauti.BeautiTesterConfig;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockType;
import dr.app.beauti.options.NucModelType;
import dr.app.beauti.options.PartitionModel;
import dr.app.beauti.options.TreePrior;


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

	public GTRParameterEstimationTest() {
		createBeastXMLFiles();
	}
       

    public void createBeastXMLFiles() {
    	BeautiOptions beautiOptions;
        BeautiTesterConfig btc = new BeautiTesterConfig();
        btc.createScriptWriter("examples/Joseph/test/run_script.sh");
    	
        for (int i = 0; i < 2; i++) {
			beautiOptions = btc.createOptions();
	        btc.importFromFile("examples/Joseph/n.nex", beautiOptions, false); 
	        
	        PartitionModel model = beautiOptions.getPartitionModels().get(0);

	        model.setNucSubstitutionModel(NucModelType.GTR);
	        model.setCodonHeteroPattern(null);
	        model.setUnlinkedSubstitutionModel(false);
	        model.setUnlinkedHeterogeneityModel(false);
	        
	        model.setGammaHetero(false);
	        model.setGammaCategories(4);
	        model.setInvarHetero(false);
	        	        
	        beautiOptions.nodeHeightPrior = TreePrior.YULE;
	        beautiOptions.clockType = ClockType.STRICT_CLOCK;
	        	
	        beautiOptions.chainLength = 2000000;
	        
	        btc.generate("examples/Joseph/test/gtr_" + Integer.toString(i) + "_", beautiOptions);        
	        
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
