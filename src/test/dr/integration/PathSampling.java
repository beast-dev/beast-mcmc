package test.dr.integration;

import dr.evolution.sequence.Sequence;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.oldevomodel.substmodel.HKY;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.math.distributions.GammaDistribution;
import dr.math.LogTricks;

/**
 * @author Marc A. Suchard
 */
public class PathSampling {

    public PathSampling() {

        int simulationLength = 10000000;
        double kappaPriorShape = 0.1;
        double kappaPriorScale = 1.0;
        double rootPriorShape = 1.0;
        double rootPriorScale = 1.0;

        // Computes the marginal likelihood of two aligned DNA sequences in a tree with height T under the Kimura model
        Sequence seq1 = new Sequence("AGAGCTCTAAGAGCTCTAAGAGCTCTAAGAGCTCTA");
        Sequence seq2 = new Sequence("AAGGCCTTTAAGGCCTTTAAGGCCTTTAAGGCCTTT");
        if (seq1.getLength() != seq2.getLength()) {
            throw new RuntimeException("Data sequences must be the same length");
        }
        DataType dataType = Nucleotides.INSTANCE;
        seq1.setDataType(dataType);
        seq2.setDataType(dataType);
        final int seqLength = seq1.getLength();
        final int stateCount = dataType.getStateCount();
        final double[] transitionProbabilities = new double[stateCount*stateCount];
        final double[] frequencies = new double[stateCount];
        for(int i = 0; i < stateCount; i++) {
            frequencies[i] = 1.0 / stateCount;
        }

        Parameter kappa = new Parameter.Default(10.0);
        HKY hky = new HKY(kappa, new FrequencyModel(dataType, frequencies));
        Parameter root = new Parameter.Default(1.0);

        GammaDistribution kappaPrior = new GammaDistribution(kappaPriorShape, kappaPriorScale);
        GammaDistribution rootPrior = new GammaDistribution(rootPriorShape, rootPriorScale);

        double logMarginalLikelihood = 0;
        for(int i = 0; i < simulationLength; i++) {

            if (i % 100000 == 0) {
                System.err.println("iter: "+i);
            }

            // Draw new parameters from prior
            kappa.setParameterValue(0,kappaPrior.nextGamma());
            root.setParameterValue(0,rootPrior.nextGamma());

            // Compute likelihood * prior
            double logPosterior = 0;
            hky.getTransitionProbabilities(root.getParameterValue(0), transitionProbabilities);
            for(int j = 0; j < seqLength; j++) {
                int char1 = seq1.getState(j);
                int char2 = seq2.getState(j);
                double prob = 0;
                for(int k = 0; k < stateCount; k++ ) {
                    prob += frequencies[k] * transitionProbabilities[k * stateCount + char1] *
                                             transitionProbabilities[k * stateCount + char2];
                }
                logPosterior += Math.log(prob);
            }

            // Add to logMarginalLikelihood on the unit-scale
            if (i == 0) {
                logMarginalLikelihood = logPosterior;
            } else {
                logMarginalLikelihood = LogTricks.logSum(logMarginalLikelihood,logPosterior);
            }
        }

        // Normalize
        logMarginalLikelihood -= Math.log(simulationLength);
        System.out.println("logMarginalLikelihood = "+logMarginalLikelihood);
    }

    public static void main(String[] args) {
        new PathSampling();
    }

}
