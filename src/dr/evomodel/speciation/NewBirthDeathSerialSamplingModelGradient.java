package dr.evomodel.speciation;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

public class NewBirthDeathSerialSamplingModelGradient implements SpeciationModelGradientProvider {

    private final NewBirthDeathSerialSamplingModel model;
    private double[] savedGradient;
    private BigFastTreeIntervals savedTreeInterval;

    private double savedQ;
    private double[] savedPartialQ;

    public NewBirthDeathSerialSamplingModelGradient(NewBirthDeathSerialSamplingModel model) {
        this.model = model;
        this.savedGradient = null;
        this.savedTreeInterval = null;
        this.savedQ = Double.MIN_VALUE;
        this.savedPartialQ = null;
    }

    // TODO(yucais): call these functions when a new tree comes!
    public void clearGradient() {
        this.savedGradient = null;
    }
    public void clearTreeInterval(){
        this.savedTreeInterval = null;
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

    private double Q(double t){
        return 4*(1/Math.exp(model.logq(t)));
    }

    // Gradient w.r.t. Rho
    private double[] partialC1C2partialRho() {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();
        double C1 = constants[0];

        double[] partialC1C2 = new double[2];
        partialC1C2[0] = 0;
        partialC1C2[1] = 2 * lambda / C1;

        return partialC1C2;
    }
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

    private double[] partialC1C2partialLambda() {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();
        double mu = model.mu();
        double psi = model.psi();
        double rho = model.rho();
        double C1 = constants[0];

        double[] partialC1C2 = new double[2];
        partialC1C2[0] = (lambda - mu + psi) / C1;
        partialC1C2[1] = ((2*rho - 1)*C1 - (-lambda + mu + 2 * lambda * rho + psi) * partialC1C2[0]) / (C1 * C1);

        return partialC1C2;
    }

    private double[] partialC1C2partialPsi() {
        // c1 == constants[0], c2 == constants[1]
        double[] constants = model.getConstants();
        double lambda = model.lambda();
        double mu = model.mu();
        double psi = model.psi();
        double rho = model.rho();
        double C1 = constants[0];

        double[] partialC1C2 = new double[2];
        partialC1C2[0] = (lambda + mu + psi) / C1;
        partialC1C2[1] = (C1 + (lambda - mu - 2 * lambda * rho - psi) * partialC1C2[0]) / (C1 * C1);

        return partialC1C2;
    }

    @Override
    public Parameter getSamplingProbabilityParameter() {
        return model.samplingFractionAtPresent;
    }

    @Override
    public double[] getSamplingProbabilityGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[3];
        return result;
    }

    @Override
    public Parameter getDeathRateParameter() {
        return model.deathRate;
    }

    @Override
    public double[] getDeathRateGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[1];
        return result;
    }


    @Override
    public Parameter getBirthRateParameter() {
        return model.birthRate;
    }

    @Override
    public double[] getBirthRateGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[0];
        return result;
    }

    @Override
    public Parameter getSamplingRateParameter() {
        return model.serialSamplingRate;
    }

    @Override
    public double[] getSamplingRateGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[2];
        return result;
    }


    @Override
    public Parameter getTreatmentProbabilityParameter() {
        return model.treatmentProbability;
    }
    @Override
    public double[] getTreatmentProbabilityGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[4];
        return result;
    }

    // gradients for all
    // (lambda, mu, psi, rho)
    public double[][] partialC1C2partialAll(){
        double[][] partialC1C2_all = new double[4][2];
        partialC1C2_all[0] = partialC1C2partialLambda();
        partialC1C2_all[1] = partialC1C2partialMu();
        partialC1C2_all[2] = partialC1C2partialPsi();
        partialC1C2_all[3] = partialC1C2partialRho();
        return partialC1C2_all;
    }

    // (lambda, mu, psi, rho)
    public double[] partialQpartialAll(double t) {
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];
        double v = Math.exp(C1 * t) * (1 + C2) - Math.exp(-C1 * t) * (1 - C2) - 2 * C2;
        double v1 = Math.exp(C1 * t) * (1 + C2) * (1 + C2) - Math.exp(-C1 * t) * (1 - C2) * (1 - C2);

        double[][] partialC1C2_all = partialC1C2partialAll();

        double[] partialQ_all = new double[4];
        for (int i = 0; i < 4; ++i) {
            partialQ_all[i] += t * partialC1C2_all[i][0] * v1;
            partialQ_all[i] += 2 * partialC1C2_all[i][1] * v;
        }
        return partialQ_all;
    }

    // (lambda, mu, psi, rho)
    public double[] partialG2partialAll(double t) {
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];

        double[][] partialC1C2_all = partialC1C2partialAll();

        double[] partialG2_all = new double[4];
        for (int i = 0; i < 3; ++i) {
            double partialC1 = partialC1C2_all[i][0];
            double partialC2 = partialC1C2_all[i][1];
            double partialG2 = g1(t) * ((partialC1 * (1 + C2) + partialC2 * C1)) - (partialC1 * t * Math.exp(-C1 * t) * (C2 - 1) + partialC2 * (1 - Math.exp(-C1 * t))) * C1 * (1 + C2);
            double G1 = g1(t);
            partialG2 = -2 * partialG2 / (G1 * G1);
            partialG2 += partialC1;
            partialG2_all[i] = partialG2;
        }
        partialG2_all[3] = 0; // w.r.t. rho
        return partialG2_all;
    }

    // (lambda, mu, psi, rho)
    public double[] partialP0partialAll(double t) {
        double[] partialG2_all = partialG2partialAll(t);

        double[] partialP0_all = new double[4];

        double lambda = model.lambda();
        double G2 = g2(t);
        double mu = model.mu();
        double psi = model.psi();
        double[] constants = model.getConstants();
        double C1 = constants[0];
        double C2 = constants[1];
        double G1 = g1(t);

        // lambda
        partialP0_all[0] = (-mu - psi + lambda * partialG2_all[0] - G2) / (2 * lambda*lambda);
        // mu
        partialP0_all[1] = (1 + partialG2_all[1]) / (2 * lambda);
        // psi
        partialP0_all[2] = (1 + partialG2_all[2]) / (2 * lambda);
        // rho
        partialP0_all[3] = -C1 / lambda * (2 * lambda / C1 * (G1 - (1 - Math.exp(-C1 * t)) * (1 + C2))) / (G1 * G1);

        return partialP0_all;
    }


    // (lambda, mu, psi, rho)
    public double[] getAllGradient(Tree tree, NodeRef node) {
        return getGradientLogDensityImpl((TreeModel) tree);
//        if (savedGradient != null) {
//            return savedGradient;
//        }
//        else {
//            model.precomputeConstants();
//
//            double[] partialLL_all = new double[4];
//
//            double origin = model.originTime.getValue(0);
//            double p0 = model.p0(origin);
//            double Q = Q(origin);
//            double lambda = model.lambda();
//            double psi = model.psi();
//            double rho = model.rho();
//
//            double[] partialP0_all_origin = partialP0partialAll(origin);
//            double[] partialQ_all_origin = partialQpartialAll(origin);
//
//
//            for (int i = 0; i < 4; ++i) {
//               // partialLL_all[i] = 1 / (1 - p0) * partialP0_all_origin[i];
//                partialLL_all[i] -= partialQ_all_origin[i] / Q;
//            }
//
//            BigFastTreeIntervals treeIntervals = this.savedTreeInterval;
//
//            if (this.savedTreeInterval == null) {
//                treeIntervals = new BigFastTreeIntervals((TreeModel) tree);
//                this.savedTreeInterval = treeIntervals;
//            }
//
//            int m = 0;
//            int mPlusn = 1;
//            int k = 0;
//
//            for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {
//                double intervalStart = treeIntervals.getIntervalTime(i);
//                final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
//                // final int nLineages = treeIntervals.getLineageCount(i);
//                if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
//                    double t = intervalEnd;
//                    Q = Q(t);
//                    double[] partialQ_all = partialQpartialAll(t);
//                    for (int j = 0; j < 4; ++j) {
//                        partialLL_all[j] -= partialQ_all[j] / Q;
//                    }
//                    mPlusn += 1;
//                } else if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {
//                    double t = intervalEnd;
//                    Q = Q(t);
//                    double[] partialP0_all = partialP0partialAll(t);
//                    double[] partialQ_all = partialQpartialAll(t);
//                    double P0 = model.p0(t);
//                    double r = model.r();
//                    double v = (1 - r) / ((1 - r) * P0 + r);
//                    for (int j = 0; j < 4; ++j) {
//                        partialLL_all[j] += v * partialP0_all[j] + partialQ_all[j] / Q;
//                    }
//                    m += 1;
//                }
//            }
//
//            // post processing
//            // (lambda, mu, psi, rho)
//            // lambda
//            int n = mPlusn - m;
//            partialLL_all[0] += (n + m - 1) / lambda;
//            // mu - skip
//            // psi
//            partialLL_all[2] += (k + m) / psi;
//            // rho
//            if (rho != 0) {
//                partialLL_all[3] += n / rho;
//            } else {
//                double[] partialP0_all = partialP0partialAll(0);
//                double P0 = model.p0(0);
//                double r = model.r();
//                double v = (1 - r) / ((1 - r) * P0 + r);
//                for (int j = 0; j < 4; ++j) {
//                    partialLL_all[j] += v * partialP0_all[j];
//                }
//                partialLL_all[2] += 1/psi;
//            }
//            savedGradient = partialLL_all;
//
//            return partialLL_all;
//        }


    }


    // @Override
    public void precomputeGradientConstants() {
        model.precomputeConstants();
        this.savedQ = Double.MIN_VALUE;
        this.savedPartialQ = null;
    }

    // @Override
    public void processGradientModelSegmentBreakPoint(double[] gradient, int currentModelSegment, double intervalStart, double segmentIntervalEnd) {
        return;
    }

    // @Override
    public void processGradientInterval(double[] gradient, int currentModelSegment, double intervalStart, double intervalEnd, int nLineages) {
        double tOld = intervalEnd;
        double tYoung = intervalStart;
        double[] partialQ_all_old = partialQpartialAll(tOld);
        double[] partialQ_all_young;
        double Q_Old = Q(tOld);
        double Q_young;
        if (this.savedQ != Double.MIN_VALUE) {
            Q_young = this.savedQ;
        }
        else {
            Q_young = Q(tYoung);
        }
        this.savedQ = Q_Old;

        if (this.savedPartialQ != null) {
            partialQ_all_young = this.savedPartialQ;
        }
        else {
            partialQ_all_young = partialQpartialAll(tYoung);
        }
        this.savedPartialQ = partialQ_all_old;

        for (int j = 0; j < 4; ++j) {
            gradient[j] += nLineages*(partialQ_all_young[j] / Q_young - partialQ_all_old[j] / Q_Old);
        }

    }

    // @Override
    public void processGradientSampling(double[] gradient, int currentModelSegment, double intervalEnd) {
        double r = model.r();
        double rho = model.rho();
        double psi = model.psi();
        double t = intervalEnd;

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = model.rho() < Double.MIN_VALUE;

        if (noSamplingAtPresent || t > timeZeroTolerance) {
            double[] partialP0_all = partialP0partialAll(t);
            double P0 = model.p0(t);
            double v = (1 - r) / ((1 - r) * P0 + r);
            for (int j = 0; j < 4; ++j) {
                gradient[j] += v * partialP0_all[j];
            }
            gradient[2] += 1 / psi;
            gradient[4] += (1 - P0) / ((1 - r)*P0 + r);
        } else {
            if (rho != 0) {
                gradient[3] += 1 / rho;
            }
        }

    }

    // @Override
    public void processGradientCoalescence(double[] gradient, int currentModelSegment, double intervalEnd) {
        gradient[0] += 1 / model.lambda();
    }

    // @Override
    public void processGradientOrigin(double[] gradient, int currentModelSegment, double totalDuration) {
        double origin = model.originTime.getValue(0);
        double[] partialQ_all_origin = partialQpartialAll(origin);
        double[] partialQ_all_root = partialQpartialAll(totalDuration);
        for (int i = 0; i < 4; ++i) {
            // partialLL_all[i] = 1 / (1 - p0) * partialP0_all_origin[i];
            gradient[i] += partialQ_all_root[i]/Q(totalDuration) - partialQ_all_origin[i] / Q(origin);
        }

    }

    // @Override
    public void logConditioningProbability(double[] gradient) {
        return;
    }

    public double[] getBreakPoints() {
        return model.getBreakPoints();
    }



    private double[] getGradientLogDensityImpl(TreeModel tree) {

        int dimension = 5; // TODO Get from delegate
        double[] gradient = new double[dimension];

        precomputeGradientConstants(); // TODO hopefully get rid of this

        // TODO Make cached class-object
        BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals(tree);

        double[] modelBreakPoints = getBreakPoints();
        assert modelBreakPoints[modelBreakPoints.length - 1] == Double.POSITIVE_INFINITY;

        int currentModelSegment = 0;


        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {

            double intervalStart = treeIntervals.getIntervalTime(i);
            final double intervalEnd = intervalStart + treeIntervals.getInterval(i);
            final int nLineages = treeIntervals.getLineageCount(i);

            while (intervalEnd > modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?

                final double segmentIntervalEnd = modelBreakPoints[currentModelSegment];
                processGradientModelSegmentBreakPoint(gradient, currentModelSegment, intervalStart, segmentIntervalEnd);
                intervalStart = segmentIntervalEnd;
                ++currentModelSegment;
            }

            // TODO Need to check for intervalStart == intervalEnd?
            // TODO Need to check for intervalStart == intervalEnd == 0.0?

            processGradientInterval(gradient, currentModelSegment, intervalStart, intervalEnd, nLineages);

            // Interval ends with a coalescent or sampling event at time intervalEnd
            if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {

                processGradientSampling(gradient, currentModelSegment, intervalEnd);

            } else if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {

                processGradientCoalescence(gradient, currentModelSegment, intervalEnd);

            } else {
                throw new RuntimeException("Birth-death tree includes non birth/death/sampling event.");
            }
        }

        // We've missed the first sample and need to add it back
        // TODO May we missed multiple samples @ t == 0.0?
        processGradientSampling(gradient, 0, treeIntervals.getStartTime()); // TODO for-loop for models with multiple segments?

        // origin branch is a fake branch that doesn't exist in the tree, now compute its contribution
        processGradientOrigin(gradient, currentModelSegment, treeIntervals.getTotalDuration());

        logConditioningProbability(gradient);

        return gradient;

    }



}
