package test.dr.app.beagle;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import test.dr.inference.trace.TraceCorrelationAssert;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.alignment.SitePatterns;
import dr.inference.model.Parameter;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.oldevomodelxml.substmodel.HKYParser;
import dr.oldevomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.math.MathUtils;

/**
 * @author Marc Suchard
 */

public class MarkovJumpsTest extends TraceCorrelationAssert {

    public MarkovJumpsTest(String name) {
        super(name);
    }

    public void testMarkovJumps() {

        MathUtils.setSeed(666);

        createAlignment(sequencesTwo, Nucleotides.INSTANCE);

        try {
//            createSpecifiedTree("((human:1,chimp:1):1,gorilla:2)");
            createSpecifiedTree("(human:1,chimp:1)");
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse Newick tree");
        }

        //substitutionModel
        Parameter freqs = new Parameter.Default(new double[]{0.40, 0.25, 0.25, 0.10});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 10.0, 0, 100);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
//        double alpha = 0.5;
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 0.5, 0, Double.POSITIVE_INFINITY);
//        Parameter pInv = new Parameter.Default("pInv", 0.5, 0, 1);
        Parameter pInv = null;
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gammaModel", mu, null, -1, pInv);
        siteRateModel.setSubstitutionModel(hky);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        BranchModel branchModel = new HomogeneousBranchModel(
                siteRateModel.getSubstitutionModel());

        BranchRateModel branchRateModel = null;

        MarkovJumpsBeagleTreeLikelihood mjTreeLikelihood = new MarkovJumpsBeagleTreeLikelihood(
                patterns,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                null,
                false,
                PartialsRescalingScheme.AUTO,
                true,
                null,
                hky.getDataType(),
                "stateTag",
                false, // use MAP
                true, // return ML
                false, // use uniformization
                false,
                1000
        );

        int nRegisters = registerValues.length;
        int nSim = 10000;

        for (int i = 0; i < nRegisters; i++) {

            Parameter registerParameter = new Parameter.Default(registerValues[i]);

            registerParameter.setId(registerTages[i]);
            mjTreeLikelihood.addRegister(
                    registerParameter,
                    registerTypes[i],
                    registerScales[i]);
        }

        double logLikelihood = mjTreeLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);

        double[] averages = new double[nRegisters];

        for (int i = 0; i < nSim; i++) {
            for (int r = 0; r < nRegisters; r++) {
                double[][] values = mjTreeLikelihood.getMarkovJumpsForRegister(treeModel, r);
                for (double[] value : values) {
                    averages[r] += value[0];
                }
            }
            mjTreeLikelihood.makeDirty();
        }

        for (int r = 0; r < nRegisters; r++) {
            averages[r]  /= (double) nSim;
            System.out.print(" " + averages[r]);
        }

        System.out.println("");
        assertEquals(valuesFromR, averages, 1E-2);
    }

    private static double[][] registerValues = {
            {
                    0, 1, 1, 1,
                    1, 0, 1, 1,
                    1, 1, 0, 1,
                    1, 1, 1, 0
            },
            {
                    0, 1, 1, 1,
                    0, 0, 1, 1,
                    0, 0, 0, 1,
                    0, 0, 0, 0
            },
            {
                    1, 0, 0, 1
            }
    };

    private static String[] registerTages = {"jump", "upper", "reward"};

    private double[] valuesFromR = { 0.782, 0.225, 1.777 }; // No mixture
//    private double[] valuesFromR = { 0.979, 0.329, 1.777 }; // 50/50 mixture of invariance

    private static MarkovJumpsType[] registerTypes =
            {MarkovJumpsType.COUNTS, MarkovJumpsType.COUNTS, MarkovJumpsType.REWARDS};

    private static boolean[] registerScales = {false, false, true};


    static private String sequencesTwo[][] = {
            {"human", "chimp"},
            {
                    "A",
                    "A"}
    };

//    static private String sequencesThree[][] = {
//            {"human", "chimp", "gorilla"},
//            {
//                    "A",
//                    "A",
//                    "N"}
//    };
}