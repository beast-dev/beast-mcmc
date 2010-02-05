package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GTR;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.HKYParser;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;


/**
 * @author Walter Xie
 * convert testLikelihood.xml in the folder /example
 */

public class LikelihoodTest extends TraceCorrelationAssert {

    private TreeModel treeModel;
    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public LikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        createAlignment(HOMINID_TAXON_SEQUENCE, Nucleotides.INSTANCE);

        createTreeModel ();
    }

    private void createTreeModel () {

        SimpleNode[] nodes = new SimpleNode[10];
        for (int n=0; n < 10; n++) {
            nodes[n] = new SimpleNode();
        }

//        nodes[0].setHeight(0);
        nodes[0].setTaxon(taxa[0]); // human

        nodes[1].setTaxon(taxa[1]); // chimp

        nodes[2].setTaxon(taxa[2]); // bonobo

        nodes[3].setHeight(0.010772);
        nodes[3].addChild(nodes[1]);
        nodes[3].addChild(nodes[2]);

        nodes[4].setHeight(0.024003);
        nodes[4].addChild(nodes[0]);
        nodes[4].addChild(nodes[3]);

        nodes[5].setTaxon(taxa[3]); // gorilla

        nodes[6].setHeight(0.036038);
        nodes[6].addChild(nodes[4]);
        nodes[6].addChild(nodes[5]);

        nodes[7].setTaxon(taxa[4]); // orangutan

        nodes[8].setHeight(0.069125);
        nodes[8].addChild(nodes[6]);
        nodes[8].addChild(nodes[7]);

        nodes[9].setTaxon(taxa[5]); // siamang

        SimpleNode root = new SimpleNode();
        root.setHeight(0.099582);
        root.addChild(nodes[8]);
        root.addChild(nodes[9]);

        Tree tree = new SimpleTree(root);
        tree.setUnits(Units.Type.YEARS);

        treeModel = new TreeModel(tree); //treeModel
    }

    public void testNewickTree() {
        System.out.println("\nTest Simple Node to convert Newick Tree:");
        String expectedNewickTree = "((((human:0.024003,(chimp:0.010772,bonobo:0.010772):0.013231):0.012035," +
                "gorilla:0.036038):0.033087,orangutan:0.069125):0.030457,siamang:0.099582);";
        
        assertEquals("Fail to covert the correct tree !!!", expectedNewickTree, Tree.Utils.newick(treeModel, 6));
    }



    public void testLikelihoodJC69() {
        System.out.println("\nTest Likelihood using JC69:");
        // Sub model
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodJC69", format.format(-1992.20564), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodK80() {
        System.out.println("\nTest Likelihood using K80:");
        // Sub model
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 27.402591, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodK80", format.format(-1856.30305), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85() {
        System.out.println("\nTest Likelihood using HKY85:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 29.739445, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);          

        assertEquals("treeLikelihoodHKY85", format.format(-1825.21317), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85G() {
        System.out.println("\nTest Likelihood using HKY85G:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 38.829740, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModel.GAMMA_SHAPE, 0.137064, 0, 1000.0);

        GammaSiteModel siteModel = new GammaSiteModel(hky, mu, shape, 4, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodHKY85G", format.format(-1789.75936), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85I() {
        System.out.println("\nTest Likelihood using HKY85I:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 38.564672, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter invar = new Parameter.Default(GammaSiteModel.PROPORTION_INVARIANT, 0.701211, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(hky, mu, null, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodHKY85I", format.format(-1789.91240), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85GI() {
        System.out.println("\nTest Likelihood using HKY85GI:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 39.464538, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModel.GAMMA_SHAPE, 0.587649, 0, 1000.0);
        Parameter invar = new Parameter.Default(GammaSiteModel.PROPORTION_INVARIANT, 0.486548, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(hky, mu, shape, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodHKY85GI", format.format(-1789.63923), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTR() {
        System.out.println("\nTest Likelihood using GTR:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTR.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTR.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTR.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTR.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTR.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTR.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, null, 4, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTR", format.format(-1969.14584), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRI() {
        System.out.println("\nTest Likelihood using GTRI:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTR.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTR.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTR.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTR.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTR.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTR.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter invar = new Parameter.Default(GammaSiteModel.PROPORTION_INVARIANT, 0.5, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, null, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTRI", format.format(-1948.84175), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRG() {
        System.out.println("\nTest Likelihood using GTRG:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTR.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTR.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTR.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTR.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTR.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTR.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModel.GAMMA_SHAPE, 0.5, 0, 100.0);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, shape, 4, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTRG", format.format(-1949.03601), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRGI() {
        System.out.println("\nTest Likelihood using GTRGI:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTR.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTR.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTR.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTR.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTR.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTR.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModel.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModel.GAMMA_SHAPE, 0.5, 0, 100.0);
        Parameter invar = new Parameter.Default(GammaSiteModel.PROPORTION_INVARIANT, 0.5, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, shape, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTRGI", format.format(-1947.58294), format.format(treeLikelihood.getLogLikelihood()));
    }


    public static Test suite() {
        return new TestSuite(LikelihoodTest.class);
    }
}
