package test.dr.math;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.evolution.datatype.Nucleotides;
import dr.inference.markovjumps.*;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
public class UniformizedStateHistoryTest extends MathTestCase {

    public void setUp() {

        MathUtils.setSeed(666);

        Parameter kappa = new Parameter.Default(1, 2.0);
        double[] pi = {0.45, 0.05, 0.30, 0.20};
        Parameter freqs = new Parameter.Default(pi);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        hky = new HKY(kappa, f);
        stateCount = hky.getDataType().getStateCount();

        double[] lambda = new double[stateCount * stateCount];
        hky.getInfinitesimalMatrix(lambda);

        process = new SubordinatedProcess(lambda, stateCount);
    }

    HKY hky;
    SubordinatedProcess process;
    int stateCount;

    public void testSubordinatedProcessGeneration() {

        double[] oneStep = process.getDtmcProbabilities(1);
        double[] threeStep = process.getDtmcProbabilities(3);

        double[] rOneStep = {
                0.2608696, 0.5217391, 4.347826e-02, 0.1739130,
                0.7826087, 0.0000000, 4.347826e-02, 0.1739130,
                0.3913043, 0.2608696, 2.220446e-16, 0.3478261,
                0.3913043, 0.2608696, 8.695652e-02, 0.2608696
        };
        MarkovJumpsCore.makeComparableToRPackage(rOneStep);

        double[] rThreeStep = {
                0.3935235, 0.3570313, 0.04988904, 0.1995562,
                0.5355470, 0.2150078, 0.04988904, 0.1995562,
                0.4490014, 0.2993343, 0.04980685, 0.2018575,
                0.4490014, 0.2993343, 0.05046437, 0.2012000
        };
        MarkovJumpsCore.makeComparableToRPackage(rThreeStep);

        System.out.println("oneStep = " + new Vector(oneStep));
        assertEquals(rOneStep, oneStep, 1E-6);

        System.out.println("threeStep = " + new Vector(threeStep));
        assertEquals(rThreeStep, threeStep, 1E-6);
    }

    public void testComputePdfForNextDraw() {

        int startState = 1;
        int endState = 0;

//        int rStartState = 3; // 1 + 1 = 2 -> 3
//        int rEndState = 1; // 0 + 1 = 1 -> 1

        double[] pdf = new double[stateCount];

        int n = 4;
        int i = 1;

        process.computePdfNextChainState(startState, endState, n, i, pdf);
        pdf = MathUtils.getNormalized(pdf);

        System.err.println("PDF = " + new Vector(pdf));

        double[] rPDF = new double[]{3.422934e-01, 3.105519e-01, 2.216160e-16, 3.471547e-01};
        MarkovJumpsCore.makeComparableToRPackage(rPDF);
        assertEquals(rPDF, pdf, 1E-6);
    }

    public void testTotalChangesSamplingMethods() {

        try {

        int startState = 1;
        int endState = 0;
        double time = 0.5;
        double[] ctmcProbabilities = new double[stateCount * stateCount];

        hky.getTransitionProbabilities(time, ctmcProbabilities);
        double ctmcProbability = ctmcProbabilities[startState * stateCount + endState];

        double[] pdf = process.computePDFDirectly(startState, endState, time, ctmcProbability, 10);

        System.out.println("PDF = " + new Vector(pdf));

        double cutoff = pdf[0] + pdf[1] - 1E-6;
        System.out.println("Test cutoff = " + cutoff);
        int draw = process.drawNumberOfChanges(startState, endState, time, ctmcProbability, cutoff);
        assertEquals(1, draw);

        cutoff = pdf[0] + pdf[1] + 1E-6;
        System.out.println("Test cutoff = " + cutoff);
        draw = process.drawNumberOfChanges(startState, endState, time, ctmcProbability, cutoff);
        assertEquals(2, draw);

        cutoff = pdf[0] + pdf[1] + pdf[2] + 1E-6;
        System.out.println("Test cutoff = " + cutoff);
        draw = process.drawNumberOfChanges(startState, endState, time, ctmcProbability, cutoff);
        assertEquals(3, draw);

        System.out.println("");

        startState = 1;
        endState = 1;
        time = 0.75;

        hky.getTransitionProbabilities(time, ctmcProbabilities);
        ctmcProbability = ctmcProbabilities[startState * stateCount + endState];
        pdf = process.computePDFDirectly(startState, endState, time, ctmcProbability, 10);

        System.out.println("PDF = " + new Vector(pdf));

        cutoff = pdf[0] + pdf[1] + pdf[2] - 1E-6;
        System.out.println("Test cutoff = " + cutoff);
        draw = process.drawNumberOfChanges(startState, endState, time, ctmcProbability, cutoff);
        assertEquals(2, draw);

        cutoff = pdf[0] + pdf[1] + pdf[2] + 1E-6;
        System.out.println("Test cutoff = " + cutoff);
        draw = process.drawNumberOfChanges(startState, endState, time, ctmcProbability, cutoff);
        assertEquals(3, draw);

        } catch (SubordinatedProcess.Exception e) {
            throw new RuntimeException("Subordinated process exception");
        }
    }

    public void testStateHistorySimulationForJumps() {

        try {

        double startingTime = 1.0;
        double endingTime = 3.0;
        int startingState = 1;
        int endingState = 3;

        int N = 1000000;

        double[] tmp = new double[stateCount * stateCount];
        hky.getTransitionProbabilities(endingTime - startingTime, tmp);
        double transitionProbability = tmp[startingState * stateCount + endingState];

        double[][] registers = new double[2][stateCount * stateCount];
        MarkovJumpsCore.fillRegistrationMatrix(registers[0], stateCount); // Count all jumps
        registers[1][2 * stateCount + 1] = 1.0; // Mark just one state!

        double[] expectations = new double[registers.length];

        for (int i = 0; i < N; i++) {

            StateHistory history = UniformizedStateHistory.simulateConditionalOnEndingState(
                    startingTime,
                    startingState,
                    endingTime,
                    endingState,
                    transitionProbability,
                    stateCount,
                    process);
            for (int j = 0; j < registers.length; j++) {
                expectations[j] += history.getTotalRegisteredCounts(registers[j]);
            }
        }

        // Determine analytic solution
        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(hky);
        double[] mjExpectations = new double[stateCount * stateCount];

        for (int j = 0; j < registers.length; j++) {
            expectations[j] /= (double) N;
            System.out.println("Expected number for register = " + expectations[j]);

            markovjumps.setRegistration(registers[j]);
            markovjumps.computeCondStatMarkovJumps(endingTime - startingTime, mjExpectations);

            assertEquals(mjExpectations[startingState * stateCount + endingState], expectations[j], 1E-2);
        }

        } catch (SubordinatedProcess.Exception e) {
            throw new RuntimeException("Subordinated process exception");
        }
    }
    public void testStateHistorySimulationForRewards() {

        try {

        double startingTime = 1.0;
        double endingTime = 3.0;
        int startingState = 1;
        int endingState = 3;

        int N = 1000000;

        double[] tmp = new double[stateCount * stateCount];
        hky.getTransitionProbabilities(endingTime - startingTime, tmp);
        double transitionProbability = tmp[startingState * stateCount + endingState];

        double[][] registers = new double[3][stateCount];
        Arrays.fill(registers[0], 1.0); // Reward all states
        registers[1][0] = 1.0; // Reward just one state!
        registers[2][3] = 1.0; // Reward just one state!

        double[] expectations = new double[registers.length];

        for (int i = 0; i < N; i++) {

            StateHistory history = UniformizedStateHistory.simulateConditionalOnEndingState(
                    startingTime,
                    startingState,
                    endingTime,
                    endingState,
                    transitionProbability,
                    stateCount,
                    process);
            for (int j = 0; j < registers.length; j++) {
                expectations[j] += history.getTotalReward(registers[j]);
            }
        }

        // Determine analytic solution
        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(hky, MarkovJumpsType.REWARDS);
        double[] mjExpectations = new double[stateCount * stateCount];

        for (int j = 0; j < registers.length; j++) {
            expectations[j] /= (double) N;
            System.out.println("Expected reward for register[" + j +"] = " + expectations[j]);

            markovjumps.setRegistration(registers[j]);
            markovjumps.computeCondStatMarkovJumps(endingTime - startingTime, mjExpectations);

            assertEquals(mjExpectations[startingState * stateCount + endingState], expectations[j], 1E-2);
        }
            
        } catch (SubordinatedProcess.Exception e) {
            throw new RuntimeException("Subordinated process exception");
        }
    }    
}

/*

# R markovjumps code

         hky = as.eigen(hky.model(2, 1, c(0.45, 0.30, 0.05, 0.20), scale = T))
         maxRate = - min(hky$rate.matrix)
         R = diag(4) + hky$rate.matrix / maxRate
         R3 = R %*% R %*% R
         R4 = R %*% R3
         (R[3,] * R3[,1]) / R4[3,1]
 */
