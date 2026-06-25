package dr.evomodel.coalescent.basta;

import dr.evolution.tree.MutableTreeModel;
import dr.inference.model.Parameter;

public class AnchoredExponentialGrowthPopulationSizeModel extends IntervalSpecificPopulationSizeModel {

    private final Parameter logPopSizesParameter;
    private final Parameter growthRateParameter;
    private final double fixedAnchorTime;
    private final MutableTreeModel treeModel;


    private final boolean inputIsLog;


    public AnchoredExponentialGrowthPopulationSizeModel(String name,
                                                        Parameter logPopSizesParameter,
                                                        Parameter growthRateParameter,
                                                        double anchorTime,
                                                        int stateCount,
                                                        int numIntervals) {
        this(name, logPopSizesParameter, growthRateParameter, anchorTime, stateCount, numIntervals, true);
    }


    public AnchoredExponentialGrowthPopulationSizeModel(String name,
                                                        Parameter logPopSizesParameter,
                                                        Parameter growthRateParameter,
                                                        double anchorTime,
                                                        int stateCount,
                                                        int numIntervals,
                                                        boolean inputIsLog) {
        super(name, logPopSizesParameter, stateCount, numIntervals);

        this.logPopSizesParameter = logPopSizesParameter;
        this.growthRateParameter = growthRateParameter;
        this.fixedAnchorTime = anchorTime;
        this.treeModel = null;
        this.inputIsLog = inputIsLog;

        addVariable(growthRateParameter);
        validateDimensions(stateCount);
    }


    public AnchoredExponentialGrowthPopulationSizeModel(String name,
                                                        Parameter logPopSizesParameter,
                                                        Parameter growthRateParameter,
                                                        MutableTreeModel treeModel,
                                                        int stateCount,
                                                        int numIntervals) {
        this(name, logPopSizesParameter, growthRateParameter, treeModel, stateCount, numIntervals, true);
    }


    public AnchoredExponentialGrowthPopulationSizeModel(String name,
                                                        Parameter logPopSizesParameter,
                                                        Parameter growthRateParameter,
                                                        MutableTreeModel treeModel,
                                                        int stateCount,
                                                        int numIntervals,
                                                        boolean inputIsLog) {
        super(name, logPopSizesParameter, stateCount, numIntervals);

        this.logPopSizesParameter = logPopSizesParameter;
        this.growthRateParameter = growthRateParameter;
        this.fixedAnchorTime = Double.NaN;
        this.treeModel = treeModel;
        this.inputIsLog = inputIsLog;

        addVariable(growthRateParameter);
        addModel(treeModel);
        validateDimensions(stateCount);
    }

    private void validateDimensions(int stateCount) {
        int etaDim = logPopSizesParameter.getDimension();
        if (!(etaDim == 1 || etaDim == stateCount)) {
            throw new IllegalArgumentException(
                    "logPopSizes dimension must be 1 or stateCount (" + stateCount +
                            "), got " + etaDim);
        }

        int rDim = growthRateParameter.getDimension();
        if (!(rDim == 1 || rDim == stateCount)) {
            throw new IllegalArgumentException(
                    "Growth rate parameter dimension must be 1 or stateCount (" + stateCount +
                            "), got " + rDim);
        }
    }

    public double getAnchorTime() {
        if (treeModel != null) {
            return treeModel.getNodeHeight(treeModel.getRoot());
        }
        return fixedAnchorTime;
    }

    public boolean isRootAnchored() {
        return treeModel != null;
    }

    private double getEta(int state) {
        int idx = logPopSizesParameter.getDimension() == 1 ? 0 : state;
        double value = logPopSizesParameter.getParameterValue(idx);
        return inputIsLog ? value : Math.log(value);
    }

    private double getR(int state) {
        int idx = growthRateParameter.getDimension() == 1 ? 0 : state;
        return growthRateParameter.getParameterValue(idx);
    }

    @Override
    protected void getBaseSizes(double[] baseSizes) {
        double tAnchor = getAnchorTime();
        for (int k = 0; k < stateCount; ++k) {
            // N_d(0) = exp(eta_d + r_d * t_anchor)
            baseSizes[k] = Math.exp(getEta(k) + getR(k) * tAnchor);
        }
    }

    @Override
    protected double calculatePopulationSizeAtTime(int state, int interval,
                                                   double time,
                                                   double intervalStartTime,
                                                   double intervalEndTime) {
        double eta = getEta(state);
        double r = getR(state);
        double tAnchor = getAnchorTime();
        // N_d(t) = exp(eta_d - r_d * (t - t_anchor))
        return Math.exp(eta - r * (time - tAnchor));
    }

    @Override
    protected double calculateIntervalIntegral(int state, int interval,
                                               double intervalStartTime,
                                               double intervalLength) {
        double eta = getEta(state);
        double r = getR(state);
        double tAnchor = getAnchorTime();
        double t1 = intervalStartTime;
        double t2 = intervalStartTime + intervalLength;

        if (Math.abs(r) > 1e-10) {
            // 1/N_d(t) = exp(-eta_d + r_d * (t - t_anchor))
            // integral from t1 to t2:
            //   (exp(-eta + r*(t2 - tAnchor)) - exp(-eta + r*(t1 - tAnchor))) / r
            double a1 = -eta + r * (t1 - tAnchor);
            double a2 = -eta + r * (t2 - tAnchor);

            double maxExp = Math.max(a1, a2);
            if (maxExp > 700.0) {
                throw new RuntimeException(String.format(
                        "Exponential overflow in anchored integral: " +
                                "r=%g, eta=%g, tAnchor=%g, t2=%g, maxExp=%g",
                        r, eta, tAnchor, t2, maxExp));
            }

            return (Math.exp(a2) - Math.exp(a1)) / r;
        } else {
            return intervalLength * Math.exp(-eta);
        }
    }

    @Override
    public PopulationSizeModelType getModelType() {
        return PopulationSizeModelType.ROOT_ANCHORED_EXPONENTIAL;
    }

    public Parameter getLogPopSizesParameter() {
        return logPopSizesParameter;
    }

    public Parameter getGrowthRateParameter() {
        return growthRateParameter;
    }

    public MutableTreeModel getTreeModel() {
        return treeModel;
    }
}
