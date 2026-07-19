package dr.evomodel.coalescent.basta;

import dr.evolution.tree.MutableTreeModel;
import dr.inference.model.Parameter;

public class SharedAncestralExponentialGrowthPopulationSizeModel extends IntervalSpecificPopulationSizeModel {

    private final Parameter logPopSizeParameter;
    private final Parameter ancestralGrowthRateParameter;
    private final Parameter growthRateParameter;
    private final double fixedAnchorTime;
    private final Parameter anchorProportionParameter;
    private final MutableTreeModel treeModel;
    private final boolean inputIsLog;
    private final boolean perDemeAnchor;

    public SharedAncestralExponentialGrowthPopulationSizeModel(String name,
                                                               Parameter logPopSizeParameter,
                                                               Parameter ancestralGrowthRateParameter,
                                                               Parameter growthRateParameter,
                                                               double anchorTime,
                                                               int stateCount,
                                                               int numIntervals,
                                                               boolean inputIsLog) {
        super(name, logPopSizeParameter, stateCount, numIntervals);

        this.logPopSizeParameter = logPopSizeParameter;
        this.ancestralGrowthRateParameter = ancestralGrowthRateParameter;
        this.growthRateParameter = growthRateParameter;
        this.fixedAnchorTime = anchorTime;
        this.anchorProportionParameter = null;
        this.treeModel = null;
        this.inputIsLog = inputIsLog;
        this.perDemeAnchor = false;

        addVariable(ancestralGrowthRateParameter);
        addVariable(growthRateParameter);
        validateDimensions(stateCount);
    }


    public SharedAncestralExponentialGrowthPopulationSizeModel(String name,
                                                               Parameter logPopSizeParameter,
                                                               Parameter ancestralGrowthRateParameter,
                                                               Parameter growthRateParameter,
                                                               Parameter anchorProportionParameter,
                                                               MutableTreeModel treeModel,
                                                               int stateCount,
                                                               int numIntervals,
                                                               boolean inputIsLog) {
        super(name, logPopSizeParameter, stateCount, numIntervals);

        this.logPopSizeParameter = logPopSizeParameter;
        this.ancestralGrowthRateParameter = ancestralGrowthRateParameter;
        this.growthRateParameter = growthRateParameter;
        this.fixedAnchorTime = Double.NaN;
        this.anchorProportionParameter = anchorProportionParameter;
        this.treeModel = treeModel;
        this.inputIsLog = inputIsLog;
        // A per-deme anchor is requested when the proportion vector has one entry per deme
        this.perDemeAnchor = anchorProportionParameter.getDimension() > 1;

        addVariable(ancestralGrowthRateParameter);
        addVariable(growthRateParameter);
        addVariable(anchorProportionParameter);
        addModel(treeModel);
        validateDimensions(stateCount);
    }

    private void validateDimensions(int stateCount) {
        if (logPopSizeParameter.getDimension() != 1) {
            throw new IllegalArgumentException(
                    "Shared anchor population size must have dimension 1 (a single shared value across demes), got "
                            + logPopSizeParameter.getDimension());
        }

        if (ancestralGrowthRateParameter.getDimension() != 1) {
            throw new IllegalArgumentException(
                    "Ancestral growth rate must have dimension 1 (a single shared r_anc), got "
                            + ancestralGrowthRateParameter.getDimension());
        }

        int rDim = growthRateParameter.getDimension();
        if (!(rDim == 1 || rDim == stateCount)) {
            throw new IllegalArgumentException(
                    "Recent (deme-specific) growth rate dimension must be 1 or stateCount (" + stateCount +
                            "), got " + rDim);
        }

        if (anchorProportionParameter != null) {
            int apDim = anchorProportionParameter.getDimension();
            if (!(apDim == 1 || apDim == stateCount)) {
                throw new IllegalArgumentException(
                        "Anchor proportion must have dimension 1 (shared anchor) or stateCount (" +
                                stateCount + ", per-deme anchor), got " + apDim);
            }
        }
    }

    public double getAnchorTime() {
        return getAnchorTime(0);
    }


    public double getAnchorTime(int state) {
        if (anchorProportionParameter != null) {
            int idx = perDemeAnchor ? state : 0;
            double proportion = anchorProportionParameter.getParameterValue(idx);
            return proportion * treeModel.getNodeHeight(treeModel.getRoot());
        }
        return fixedAnchorTime;
    }

    public boolean isPerDemeAnchor() {
        return perDemeAnchor;
    }

    public boolean isEstimatedAnchor() {
        return anchorProportionParameter != null;
    }

    private double getEta() {
        double value = logPopSizeParameter.getParameterValue(0);
        return inputIsLog ? value : Math.log(value);
    }

    private double getRecentR(int state) {
        int idx = growthRateParameter.getDimension() == 1 ? 0 : state;
        return growthRateParameter.getParameterValue(idx);
    }

    private double getAncestralR() {
        return ancestralGrowthRateParameter.getParameterValue(0);
    }

    @Override
    protected void getBaseSizes(double[] baseSizes) {
        double eta = getEta();
        if (!perDemeAnchor) {
            double tAnchor = getAnchorTime(0);
            for (int k = 0; k < stateCount; ++k) {
                // present (t = 0) is in the recent epoch: N_d(0) = exp(eta + r_d * t_anchor)
                baseSizes[k] = Math.exp(eta + getRecentR(k) * tAnchor);
            }
        } else {
            double rAnc = getAncestralR();
            for (int k = 0; k < stateCount; ++k) {
                double tA = getAnchorTime(k);
                baseSizes[k] = Math.exp(eta + (getRecentR(k) - rAnc) * tA);
            }
        }
    }

    @Override
    protected double calculatePopulationSizeAtTime(int state, int interval,
                                                   double time,
                                                   double intervalStartTime,
                                                   double intervalEndTime) {
        double eta = getEta();
        if (!perDemeAnchor) {
            double tAnchor = getAnchorTime(0);
            double r = (time <= tAnchor) ? getRecentR(state) : getAncestralR();
            return Math.exp(eta - r * (time - tAnchor));
        }

        double tA = getAnchorTime(state);
        double rAnc = getAncestralR();
        if (time <= tA) {
            return Math.exp(eta - rAnc * tA - getRecentR(state) * (time - tA));
        }

        return Math.exp(eta - rAnc * time);
    }

    @Override
    protected double calculateIntervalIntegral(int state, int interval,
                                               double intervalStartTime,
                                               double intervalLength) {
        double eta = getEta();
        double t1 = intervalStartTime;
        double t2 = intervalStartTime + intervalLength;

        double rRecent = getRecentR(state);
        double rAncestral = getAncestralR();

        if (!perDemeAnchor) {
            double tAnchor = getAnchorTime(0);
            if (t1 < tAnchor && tAnchor < t2) {
                return segmentIntegral(eta, rRecent, tAnchor, t1, tAnchor)
                        + segmentIntegral(eta, rAncestral, tAnchor, tAnchor, t2);
            } else if (t2 <= tAnchor) {
                return segmentIntegral(eta, rRecent, tAnchor, t1, t2);
            } else {
                return segmentIntegral(eta, rAncestral, tAnchor, t1, t2);
            }
        }


        double tA = getAnchorTime(state);
        double recentC0 = -eta + (rAncestral - rRecent) * tA;
        if (t1 < tA && tA < t2) {
            return intExp(recentC0, rRecent, t1, tA)
                    + intExp(-eta, rAncestral, tA, t2);
        } else if (t2 <= tA) {
            return intExp(recentC0, rRecent, t1, t2);
        } else {
            return intExp(-eta, rAncestral, t1, t2);
        }
    }


    private double intExp(double c0, double c1, double lo, double hi) {
        if (Math.abs(c1) > 1e-10) {
            double a1 = c0 + c1 * lo;
            double a2 = c0 + c1 * hi;
            double maxExp = Math.max(a1, a2);
            if (maxExp > 700.0) {
                throw new RuntimeException(String.format(
                        "Exponential overflow in per-deme shared-ancestral integral: " +
                                "c0=%g, c1=%g, lo=%g, hi=%g, maxExp=%g", c0, c1, lo, hi, maxExp));
            }
            return (Math.exp(a2) - Math.exp(a1)) / c1;
        }
        return (hi - lo) * Math.exp(c0);
    }


    private double segmentIntegral(double eta, double r, double tAnchor, double lo, double hi) {
        if (Math.abs(r) > 1e-10) {
            double a1 = -eta + r * (lo - tAnchor);
            double a2 = -eta + r * (hi - tAnchor);

            double maxExp = Math.max(a1, a2);
            if (maxExp > 700.0) {
                throw new RuntimeException(String.format(
                        "Exponential overflow in shared-ancestral integral: " +
                                "r=%g, eta=%g, tAnchor=%g, lo=%g, hi=%g, maxExp=%g",
                        r, eta, tAnchor, lo, hi, maxExp));
            }

            return (Math.exp(a2) - Math.exp(a1)) / r;
        } else {
            return (hi - lo) * Math.exp(-eta);
        }
    }

    @Override
    public PopulationSizeModelType getModelType() {
        return PopulationSizeModelType.SHARED_ANCESTRAL_EXPONENTIAL;
    }

    public Parameter getLogPopSizeParameter() {
        return logPopSizeParameter;
    }

    public Parameter getAncestralGrowthRateParameter() {
        return ancestralGrowthRateParameter;
    }

    public Parameter getGrowthRateParameter() {
        return growthRateParameter;
    }

    public Parameter getAnchorProportionParameter() {
        return anchorProportionParameter;
    }

    public MutableTreeModel getTreeModel() {
        return treeModel;
    }
}
