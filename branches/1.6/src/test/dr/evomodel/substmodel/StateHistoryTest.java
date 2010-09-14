package test.dr.evomodel.substmodel;

import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.datatype.Nucleotides;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.StateHistory;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */
public class StateHistoryTest extends MathTestCase {

    public static final int N = 1000000;

    public void setUp() {

        MathUtils.setSeed(666);

        freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.45, 0.25, 0.05, 0.25});
        baseModel = new HKY(2.0, freqModel);
        stateCount = baseModel.getDataType().getStateCount();

        lambda = new double[stateCount * stateCount];
        baseModel.getInfinitesimalMatrix(lambda);
        System.out.println("lambda = " + new Vector(lambda));

        markovjumps = new MarkovJumpsSubstitutionModel(baseModel);
    }

    public void testFreqDistribution() {

        System.out.println("Start of FreqDistribution test");
        int startingState = 0;
        double duration = 10; // 10 expected substitutions is close to \infty

        double[] freq = new double[stateCount];

        for (int i = 0; i < N; i++) {
            StateHistory simultant = StateHistory.simulateUnconditionalOnEndingState(0.0, startingState, duration,
                    lambda, stateCount);
            freq[simultant.getEndingState()]++;
        }

        for (int i = 0; i < stateCount; i++) {
            freq[i] /= N;
        }

        System.out.println("freq = " + new Vector(freq));
        assertEquals(freq, freqModel.getFrequencies(), 1E-3);
        System.out.println("End of FreqDistribution test\n");
    }

    public void testCounts() {

        System.out.println("State of Counts test");
        int startingState = 2;
        double duration = 0.5;

        int[] counts = new int[stateCount * stateCount];
        double[] expectedCounts = new double[stateCount * stateCount];

        for (int i = 0; i < N; i++) {
            StateHistory simultant = StateHistory.simulateUnconditionalOnEndingState(0.0, startingState, duration,
                    lambda, stateCount);
            simultant.accumulateSufficientStatistics(counts, null);
        }

        for (int i = 0; i < stateCount * stateCount; i++) {
            expectedCounts[i] = (double) counts[i] / (double) N;
        }

        double[] r = new double[stateCount * stateCount];
        double[] joint = new double[stateCount * stateCount];
        double[] analytic = new double[stateCount * stateCount];

        for (int from = 0; from < stateCount; from++) {
            for (int to = 0; to < stateCount; to++) {
                double marginal = 0;
                if (from != to) {
                    MarkovJumpsCore.fillRegistrationMatrix(r, from, to, stateCount);
                    markovjumps.setRegistration(r);
                    markovjumps.computeJointStatMarkovJumps(duration, joint);

                    for (int j = 0; j < stateCount; j++) {
                        marginal += joint[startingState * stateCount + j]; // Marginalize out ending state
                    }
                }
                analytic[from * stateCount + to] = marginal;
            }
        }

        System.out.println("unconditional expected counts = " + new Vector(expectedCounts));
        System.out.println("analytic               counts = " + new Vector(analytic));

        assertEquals(expectedCounts, analytic, 1E-3);
        System.out.println("End of Counts test\n");
    }
    
    double[] lambda;
    FrequencyModel freqModel;
    SubstitutionModel baseModel;
    MarkovJumpsSubstitutionModel markovjumps;
    int stateCount;


}
