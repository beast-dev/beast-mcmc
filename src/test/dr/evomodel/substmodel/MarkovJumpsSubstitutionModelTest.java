package test.dr.evomodel.substmodel;

import test.dr.math.MathTestCase;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.evolution.datatype.Nucleotides;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */

public class MarkovJumpsSubstitutionModelTest extends MathTestCase {

    public void testMarkovJumpsCounts() {
        HKY substModel = new HKY(2.0,
                new FrequencyModel(Nucleotides.INSTANCE,
                        new double[]{0.3, 0.2, 0.25, 0.25})); // A,C,G,T

        int states = substModel.getDataType().getStateCount();

        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(substModel,
                MarkovJumpsType.COUNTS);
        double[] r = new double[states * states];
        double[] q = new double[states * states];
        double[] j = new double[states * states];
        double[] c = new double[states * states];
        double[] p = new double[states * states];

        double time = 1.0;
        int from = 0; // A
        int to = 1; // C
        MarkovJumpsCore.fillRegistrationMatrix(r, from, to, states, 1.0);
        markovjumps.setRegistration(r);

        substModel.getInfinitesimalMatrix(q);

        substModel.getTransitionProbabilities(time, p);

        markovjumps.computeJointStatMarkovJumps(time, j);

        markovjumps.computeCondStatMarkovJumps(time, c);

        MarkovJumpsCore.makeComparableToRPackage(q);
        System.out.println("Q = " + new Vector(q));

        MarkovJumpsCore.makeComparableToRPackage(p);
        System.out.println("P = " + new Vector(p));

        System.out.println("Counts:");
        MarkovJumpsCore.makeComparableToRPackage(r);
        System.out.println("R = " + new Vector(r));

        MarkovJumpsCore.makeComparableToRPackage(j);
        System.out.println("J = " + new Vector(j));

        assertEquals(rMarkovJumpsJ, j, tolerance);

        MarkovJumpsCore.makeComparableToRPackage(c);
        System.err.println("C = " + new Vector(c));

        assertEquals(rMarkovJumpsC, c, tolerance);
    }

    public void testMarkovJumpsReward() {
        HKY substModel = new HKY(2.0,
                new FrequencyModel(Nucleotides.INSTANCE,
                        new double[]{0.3, 0.2, 0.25, 0.25})); // A,C,G,T

        int states = substModel.getDataType().getStateCount();

        double[] j = new double[states * states];
        double[] c = new double[states * states];

        double time = 1.0;

        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(substModel, MarkovJumpsType.REWARDS);
        double[] rewards = {1.0, 1.0, 1.0, 1.0};
        markovjumps.setRegistration(rewards);

        markovjumps.computeJointStatMarkovJumps(time, j);

        markovjumps.computeCondStatMarkovJumps(time, c);

        System.out.println("Rewards:");
        MarkovJumpsCore.makeComparableToRPackage(rewards);
        System.out.println("R = " + new Vector(rewards));

        MarkovJumpsCore.makeComparableToRPackage(j);
        System.out.println("J = " + new Vector(j));

        assertEquals(rMarkovRewardsJ, j, tolerance);

        MarkovJumpsCore.makeComparableToRPackage(c);
        System.out.println("C = " + new Vector(c));

        assertEquals(rMarkovRewardsC, c, tolerance);
    }

    public void testMarginalRates() {
         HKY substModel = new HKY(2.0,
                new FrequencyModel(Nucleotides.INSTANCE,
                        new double[]{0.3, 0.2, 0.25, 0.25})); // A,C,G,T

        int states = substModel.getDataType().getStateCount();

        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(substModel,
                MarkovJumpsType.COUNTS);
        double[] r = new double[states * states];

        int from = 0; // A
        int to = 1; // C
        MarkovJumpsCore.fillRegistrationMatrix(r, from, to, states, 1.0);
        markovjumps.setRegistration(r);

        double marginalRate = markovjumps.getMarginalRate();
        System.out.println("Marginal rate = " + marginalRate);

        assertEquals(rMarkovMarginalRate, marginalRate, tolerance);

        MarkovJumpsCore.fillRegistrationMatrix(r, states);
        markovjumps.setRegistration(r);

        marginalRate = markovjumps.getMarginalRate();
        System.out.println("Marginal rate = " + marginalRate);
        assertEquals(1.0, marginalRate, tolerance);
    }

    private static double tolerance = 1E-6;

    private static double[] rMarkovJumpsJ = {
            0.016780099, 0.013983416, 0.08448885, 0.022470093,
            0.003179188, 0.002649323, 0.02556076, 0.004475270,
            0.001889475, 0.001574562, 0.01612631, 0.002673292,
            0.001889475, 0.001574562, 0.01612631, 0.002673292
    };

    private static double[] rMarkovJumpsC = {
            0.034557323, 0.061024826, 0.66635307, 0.141775072,
            0.011561873, 0.006024690, 0.20159458, 0.028236719,
            0.009934702, 0.009934702, 0.03850175, 0.011499342,
            0.009934702, 0.009934702, 0.08671048, 0.005744804
    };

    private static double[] rMarkovRewardsJ = {
            0.4855729, 0.2291431, 0.1267929, 0.1584911,
            0.2749717, 0.4397443, 0.1267929, 0.1584911,
            0.1901894, 0.1584911, 0.4188461, 0.2324734,
            0.1901894, 0.1584911, 0.1859787, 0.4653407
    };

    private static double[] rMarkovRewardsC = {
            1, 1, 1, 1,
            1, 1, 1, 1,
            1, 1, 1, 1,
            1, 1, 1, 1
    };

    private static double rMarkovMarginalRate = 0.2010050 * 0.3;
}

/*
# R script to compare main to original package:

subst.model.eigen = as.eigen.hky(c(2,1),
	c(0.3,0.25,0.2,0.25),scale=T)

time = 1.0
from = 0 # A
to   = 2 # C
R = matrix(0,nrow=4,ncol=4)
R[from + to*4 + 1] = 1.0

P = matexp.eigen(subst.model.eigen,time)
J = joint.mean.markov.jumps(subst.model.eigen,R,time)
C = cond.mean.markov.jumps(subst.model.eigen,R,time)

rR = c(1,1,1,1)

rJ = joint.mean.markov.rewards(subst.model.eigen,rR,time)

*/

