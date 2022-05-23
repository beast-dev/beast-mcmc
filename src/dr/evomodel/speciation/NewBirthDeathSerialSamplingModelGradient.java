package dr.evomodel.speciation;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

public class NewBirthDeathSerialSamplingModelGradient implements SpeciationModelGradientProvider {

    private final NewBirthDeathSerialSamplingModel model;

    public NewBirthDeathSerialSamplingModelGradient(NewBirthDeathSerialSamplingModel model) {
        this.model = model;
    }

    private double g1(double t) {
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];
        double G1 = Math.exp(-C1 * t) * (1 - C2) + (1 + C2);
        return G1;
    }

    private double g2(double t) {
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];
        double G1 = g1(t);
        double G2 = C1 * (1 - 2 * (1 + C2) / G1);
        return G2;
    }


    // Gradient w.r.t. Rho
    private double partialQpartialRho(double t) {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();
        double C1 = constants[0];
        double C2 = constants[1];

        // d/dc2
        double partialC2 = 2 * lambda / constants[0];

        // double partialQ = -4.0 * constants[1] * partialC2 + Math.exp(constants[0] * t )
        //         * -2.0 * (1 - constants[1]) * partialC2
        //         + Math.exp(constants[0] * t) * 2.0 * (1 + constants[1]) * partialC2;

        double partialQ = 2 * partialC2 * (Math.exp(C1 * t) * (1 + C2) - Math.exp(-C1 * t) * (1 - C2) - 2 * C2);

        return partialQ;
    }

    private double partialP0partialRho(double t) {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();
        double C1 = constants[0];
        double C2 = constants[1];

        // d/dc2
        double partialC2 = 2 * lambda / constants[0];

        // define g1
        double G1 = g1(t);
        double partialP0 = -C1 / lambda * (partialC2 * (G1 - (1 - Math.exp(-C1 * t)) * (1 + C2))) / (G1 * G1);

        return partialP0;
    }

    @Override
    public Parameter getSamplingProbabilityParameter() {
        return model.samplingFractionAtPresent;
    }

    @Override
    public double[] getSamplingProbabilityGradient(Tree tree, NodeRef node) {
        model.precomputeConstants();

        // c1 == constants[0], c2 == constants[1]
        double rho = model.rho();

        // d/dq
        double partialLL = 0.0;

        double origin = model.originTime.getValue(0);
        //double p0 = model.p0(origin);
       // partialLL += 1 / (1 - p0) * partialP0partialRho(origin);
        double Q = Math.exp(model.logq(origin));
        partialLL -= partialQpartialRho(origin) / Q;

        BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals((TreeModel) tree);

        int m = 0;
        int mPlusn = 1;


        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {
            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            // final int nLineages = treeIntervals.getLineageCount(i);

            if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {
                double t = intervalEnd;
                Q = Math.exp(model.logq(t));
                double P0 = model.p0(t);
                double r = model.r();
                partialLL += (1 - r) / ((1 - r) * P0 + r) * partialP0partialRho(t) + partialQpartialRho(t) / Q;
                m += 1;

            } else if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                double t = intervalEnd;
                Q = Math.exp(model.logq(t));
                partialLL -= partialQpartialRho(t) / Q;
                mPlusn += 1;
            }
        }

        int n = mPlusn - m;
        if (rho != 0) {
            partialLL += n / rho;
        }

        double[] partialLL_vec = new double[1];
        partialLL_vec[0] = partialLL;

        return partialLL_vec;
    }

    // Gradient w.r.t. mu
    private double[] partialC1C2partialMu() {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();
        double mu = model.mu();
        double psi = model.psi();
        double rho = model.rho();
        double C1 = constants[0];

        double[] partialC1C2 = new double[2];
        partialC1C2[0] = (-lambda + mu + psi) / C1;
        partialC1C2[1] = (C1 + (lambda - mu - 2 * lambda * rho - psi) * partialC1C2[0]) / (C1 * C1);

        return partialC1C2;
    }

    private double partialQpartialMu(double t) {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];

        double[] partialC1C2 = partialC1C2partialMu();
        double partialC1 = partialC1C2[0];
        double partialC2 = partialC1C2[1];
        double partialQ = t * partialC1 * ((Math.exp(C1 * t) * (1 + C2) * (1 + C2) - Math.exp(-C1 * t) * (1 - C2) * (1 - C2)));
        partialQ += 2 * partialC2 * ((Math.exp(C1 * t) * (1 + C2) - Math.exp(-C1 * t) * (1 - C2) - 2 * C2));
        return partialQ;
    }


    private double partialG2partialMu(double t) {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];

        double[] partialC1C2 = partialC1C2partialMu();
        double partialC1 = partialC1C2[0];
        double partialC2 = partialC1C2[1];

        double partialG2 = g1(t) * ((partialC1 * (1 + C2) + partialC2 * C1)) - (partialC1 * t * Math.exp(-C1 * t) * (C2 - 1) + partialC2 * (1 - Math.exp(-C1 * t))) * C1 * (1 + C2);
        double G1 = g1(t);
        partialG2 = -2 * partialG2 / (G1 * G1);
        partialG2 += partialC1;
        return partialG2;
    }

    private double partialP0partialMu(double t) {
        double lambda = model.lambda();
        double partialG2 = partialG2partialMu(t);
        return (1 + partialG2) / (2 * lambda);
    }

    @Override
    public Parameter getDeathRateParameter() {
        return model.deathRate;
    }

    @Override
    public double[] getDeathRateGradient(Tree tree, NodeRef node) {
        model.precomputeConstants();
        double partialLL = 0.0;
        double origin = model.originTime.getValue(0);
        // double p0 = model.p0(origin);
        //partialLL += 1 / (1 - p0) * partialP0partialMu(origin);
        double Q = Math.exp(model.logq(origin));
        partialLL -= partialQpartialMu(origin) / Q;
        BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals((TreeModel) tree);


        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {
            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            // final int nLineages = treeIntervals.getLineageCount(i);

            if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                double t = intervalEnd;
                Q = Math.exp(model.logq(t)); // can store q's when calculating likelihood
                partialLL -= partialQpartialMu(t) / Q;
            } else if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {
                double t = intervalEnd;
                Q = Math.exp(model.logq(t));
                double P0 = model.p0(t);
                double r = model.r();
                partialLL += (1 - r) / ((1 - r) * P0 + r) * partialP0partialMu(t) + partialQpartialMu(t) / Q;
            }
        }

        double[] partialLL_vec = new double[1];
        partialLL_vec[0] = partialLL;

        return partialLL_vec;
    }

}
