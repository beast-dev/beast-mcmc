package test.dr.evomodel.operators;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.evomodel.tree.TreeHeightStatistic;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreelengthStatistic;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inference.prior.Prior;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Alexei Drummond
 */
public class ExchangeOperatorTest extends TestCase {

    static final String TL = "TL";
    static final String TREE_HEIGHT = "rootHeight";

    private FlexibleTree tree5;
    private FlexibleTree tree6;

    public static Test suite() {
        return new TestSuite(ExchangeOperatorTest.class);
    }

    public void setUp() throws Exception {
        super.setUp();

        NewickImporter importer = new NewickImporter(
                "((((A:1.0,B:1.0):1.0,C:2.0),D:3.0):1.0, E:4.0);");
        tree5 = (FlexibleTree) importer.importTree(null);

        importer = new NewickImporter(
                "(((((A:1.0,B:1.0):1.0,C:2.0),D:3.0):1.0, E:4.0),F:5.0);");
        tree6 = (FlexibleTree) importer.importTree(null);
    }

    // 5 taxa trees should sample all 105 topologies
    public void testIrreducibility5() throws IOException, Importer.ImportException {
        irreducibilityTester(tree5, 105, 200000, 10);
    }

    // 6 taxa trees should sample all 945 topologies
    public void testIrreducibility6() throws IOException, Importer.ImportException {
        irreducibilityTester(tree6, 945, 2000000, 4);
    }

    public void testWideExchangeOperator2() throws IOException, ImportException {

        // probability of picking (A,B) node is 1/(2n-2) = 1/8
        // probability of swapping with D is 1/2
        // total = 1/16

        //probability of picking {D} node is 1/(2n-2) = 1/8
        //probability of picking {A,B} is 1/5
        // total = 1/40

        //total = 1/16 + 1/40 = 0.0625 + 0.025 = 0.0875
    	
    	System.out.println("Test 1: Forward");

        String treeMatch = "(((D,C),(A,B)),E);";
        
        int count = 0;
        int reps = 1000000;

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree5);
                ExchangeOperator operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
                operator.doOperation();

                String tree = Tree.Utils.newickNoLengths(treeModel);

                if (tree.equals(treeMatch)) {
                    count += 1;
                }

            } catch (OperatorFailedException e) {
                e.printStackTrace();
            }

        }
        double p_1 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_1);
        System.out.println("Number of expected ratio:\t" + 0.0875);
        assertExpectation(0.0875, p_1, reps);
        
        // since this operator is supposed to be symmetric it got a hastings ratio of one
        // this means, it should propose the same move just backwards with the same probability
        
        // BUT:
        
        // (((D:2.0,C:2.0):1.0,(A:1.0,B:1.0):2.0):1.0,E:4.0) -> ((((A,B),C),D),E)
        
        // probability of picking (A,B) node is 1/(2n-2) = 1/8
        // probability of swapping with D is 1/3
        // total = 1/24

        //probability of picking {D} node is 1/(2n-2) = 1/8
        //probability of picking {A,B} is 1/4
        // total = 1/32

        //total = 1/24 + 1/32 = 7/96 = 0.07291666666
        
    	System.out.println("Test 2: Backward");
        
        treeMatch = "((((A,B),C),D),E);";
        NewickImporter importer = new NewickImporter("(((D:2.0,C:2.0):1.0,(A:1.0,B:1.0):2.0):1.0,E:4.0);");
        FlexibleTree tree5_2 = (FlexibleTree) importer.importTree(null);

        count = 0;

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree5_2);
                ExchangeOperator operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
                operator.doOperation();

                String tree = Tree.Utils.newickNoLengths(treeModel);

                if (tree.equals(treeMatch)) {
                    count += 1;
                }

            } catch (OperatorFailedException e) {
                e.printStackTrace();
            }

        }
        double p_2 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_2);
        System.out.println("Number of expected ratio:\t" + 0.0791666);
        assertExpectation(0.0791666, p_2, reps);
    }

    /**
     * @param ep    the expected (binomial) probability of success
     * @param ap    the actual proportion of successes
     * @param count the number of attempts
     */
    protected void assertExpectation(double ep, double ap, int count) {

        if (count * ap < 5 || count * (1 - ap) < 5) throw new IllegalArgumentException();

        double stdev = Math.sqrt(ap * (1.0 - ap) * count) / count;
        double upper = ap + 2 * stdev;
        double lower = ap - 2 * stdev;

        assertTrue("Expected p=" + ep + " but got " + ap + " +/- " + stdev,
                upper > ep && lower < ep);

    }

    private void irreducibilityTester(Tree tree, int numTopologies, int chainLength, int sampleTreeEvery)
            throws IOException, Importer.ImportException {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(chainLength);
        options.setUseCoercion(true);
        options.setPreBurnin(100);
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        TreeModel treeModel = new TreeModel("treeModel", tree);
        TreelengthStatistic tls = new TreelengthStatistic(TL, treeModel);
        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);

        OperatorSchedule schedule = getWideExchangeSchedule(treeModel);

        Parameter b = new Parameter.Default("b", 2.0, 0.0, Double.MAX_VALUE);
        Parameter d = new Parameter.Default("d", 0.0, 0.0, Double.MAX_VALUE);

        SpeciationModel speciationModel = new BirthDeathGernhard08Model(b, d, Units.Type.YEARS);
        Likelihood likelihood = new SpeciationLikelihood(treeModel, speciationModel, "yule.like");

        MCLogger[] loggers = new MCLogger[2];
//        loggers[0] = new MCLogger(new ArrayLogFormatter(false), 100, false);
//        loggers[0].add(likelihood);
//        loggers[0].add(rootHeight);
//        loggers[0].add(tls);

        loggers[0] = new MCLogger(new TabDelimitedFormatter(System.out), 10000, false);
        loggers[0].add(likelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(tls);

        File file = new File("yule.trees");
        file.deleteOnExit();
        FileOutputStream out = new FileOutputStream(file);

        loggers[1] = new TreeLogger(treeModel, new TabDelimitedFormatter(out), sampleTreeEvery, true, true, false);

        mcmc.setShowOperatorAnalysis(true);

        mcmc.init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers);

        mcmc.run();
        out.flush();
        out.close();

        Set<String> uniqueTrees = new HashSet<String>();

        NexusImporter importer = new NexusImporter(new FileReader(file));
        while (importer.hasTree()) {
            Tree t = importer.importNextTree();
            uniqueTrees.add(Tree.Utils.uniqueNewick(t, t.getRoot()));
        }

        TestCase.assertEquals(numTopologies, uniqueTrees.size());
    }

    // STATIC METHODS

    public static OperatorSchedule getWideExchangeSchedule(TreeModel treeModel) {

        Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);

        ExchangeOperator operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
        ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, CoercionMode.COERCION_ON, 1.0);
        UniformOperator uniformOperator = new UniformOperator(internalHeights, 1.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(operator);
        schedule.addOperator(scaleOperator);
        schedule.addOperator(uniformOperator);

        return schedule;
    }
}
