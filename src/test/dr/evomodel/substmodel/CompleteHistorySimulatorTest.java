package test.dr.evomodel.substmodel;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.CodonLabeling;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.app.beagle.tools.CompleteHistorySimulator;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */

public class CompleteHistorySimulatorTest extends MathTestCase {

    public int N = 1000;

    public void setUp() throws Exception {
        super.setUp();
        MathUtils.setSeed(666);

        NewickImporter importer = new NewickImporter("(1:2.0,(2:1.0,3:1.0):1.0);");

        tree = importer.importTree(null);
        treeModel = new TreeModel("treeModel", tree);
    }

    public void testHKYSimulation() {

        Parameter kappa = new Parameter.Default(1, 2.0);
        double[] pi = {0.45, 0.05, 0.25, 0.25};
        Parameter freqs = new Parameter.Default(pi);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);
        int stateCount = hky.getDataType().getStateCount();

        Parameter mu = new Parameter.Default(1, 0.5);
        Parameter alpha = new Parameter.Default(1, 0.5);
        GammaSiteRateModel siteModel = new GammaSiteRateModel("gammaModel", mu, alpha, 4, null);
        siteModel.setSubstitutionModel(hky);
        BranchRateModel branchRateModel = new DefaultBranchRateModel();

        double analyticResult = TreeUtils.getTreeLength(tree, tree.getRoot()) * mu.getParameterValue(0);
        int nSites = 200;

        double[] register1 = new double[stateCount * stateCount];
        double[] register2 = new double[stateCount * stateCount];

        MarkovJumpsCore.fillRegistrationMatrix(register1, stateCount); // Count all jumps

        // Move some jumps from 1 to 2
        register1[1 * stateCount + 2] = 0;
        register2[1 * stateCount + 2] = 1;

        register1[1 * stateCount + 3] = 0;
        register2[1 * stateCount + 3] = 1;

        register1[2 * stateCount + 3] = 0;
        register2[2 * stateCount + 3] = 1;

        runSimulation(tree, siteModel, branchRateModel, nSites,
                new double[][] {register1, register2}, analyticResult);
    }

    public void testCodonSimulation() {
        Parameter kappa = new Parameter.Default(1, 2.0);
        Parameter omega = new Parameter.Default(1, 5.0); // Expect many more non-syn changes

        Codons codons = Codons.UNIVERSAL;
        int stateCount = codons.getStateCount();
        double[] p = new double[stateCount];
        for (int i = 0; i < stateCount; i++) {
            p[i] = 1.0 / (double) stateCount;
        }
        Parameter pi = new Parameter.Default(p);
        FrequencyModel f = new FrequencyModel(codons, pi);
        GY94CodonModel codonModel = new GY94CodonModel(codons, omega, kappa, f);

        Parameter mu = new Parameter.Default(1, 0.5);
        Parameter alpha = new Parameter.Default(1, 0.5);
        GammaSiteRateModel siteModel = new GammaSiteRateModel("gammaModel", mu, alpha, 4, null);
        siteModel.setSubstitutionModel(codonModel);
        BranchRateModel branchRateModel = new DefaultBranchRateModel();

        double analyticResult = TreeUtils.getTreeLength(tree, tree.getRoot()) * mu.getParameterValue(0);
        int nSites = 100;

        double[] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, false); // use base 61
        double[] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, false); // use base 61

        runSimulation(tree, siteModel, branchRateModel, nSites,
                new double[][] {synRegMatrix, nonSynRegMatrix}, analyticResult);

    }

    protected void runSimulation(Tree tree, GammaSiteRateModel siteModel, BranchRateModel branchRateModel, int nSites,
                               double[][] registers, double analyticResult) {
        runSimulation(N, tree, siteModel, branchRateModel, nSites, registers, analyticResult, null, null);
    }

    protected void runSimulation(int N, Tree tree, GammaSiteRateModel siteModel, BranchRateModel branchRateModel, int nSites,
                               double[][] registers, double analyticResult, Parameter variableParam,
                               Parameter valuesParam) {

        CompleteHistorySimulator simulator = new CompleteHistorySimulator(tree, siteModel, branchRateModel, nSites,
                false, variableParam, valuesParam);

        for (int r = 0; r < registers.length; r++) {
            Parameter registerParameter = new Parameter.Default(registers[r]);
            registerParameter.setId("set" + (r+1));
            simulator.addRegister(registerParameter, MarkovJumpsType.COUNTS, false);
        }

        int nJumps = simulator.getNumberOfJumpProcess();
        double[] results = new double[nJumps];

        for (int rep = 0; rep < N; rep++) {
            simulator.simulate(); // Generate new history

            for (int jump = 0; jump < nJumps; jump++) {
                results[jump] += accumulateJumpsOverTree(tree, simulator, jump);
            }
        }

        for (int jump = 0; jump < nJumps; jump++) {
            results[jump] /= (double) N * (double) nSites;
        }

        if (nJumps > 0) {
            System.out.println("simulations = " +  new Vector(results));
        } else {
            System.out.println("No jump processes!");
        }

        double simulationResult = accumulate(results);

        assertEquals(analyticResult, simulationResult, 1E-2);
    }

    private double accumulateJumpsOverTree(Tree tree, CompleteHistorySimulator sim, int jump) {
        double total = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) { // Sum over tree
            double[] realizedCounts = sim.getMarkovJumpsForNodeAndRegister(tree, tree.getNode(i), jump);
            total += accumulate(realizedCounts); // Sum over sites
        }
        return total;
    }

    Tree tree;
    TreeModel treeModel;

}
