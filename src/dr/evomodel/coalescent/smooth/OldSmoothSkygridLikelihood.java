package dr.evomodel.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.TreeIntervals;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.NumericalDerivative;
import dr.math.UnivariateFunction;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

import java.util.Arrays;
import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class OldSmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Citable {

    private final List<IntervalList> intervalList;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;
    private final double[] gridPoints;
    private final int numGridIntervals;



    private Type units;

    private double[] coalescentCount;
    private double[] coalescentIntensity;

    private double[] savedCoalescentCount;
    private double[] savedCoalescentIntensity;

    private final static UnivariateRealIntegrator integrator = new RombergIntegrator();
//    private final static UnivariateRealIntegrator integrator = new SimpsonIntegrator();

    public OldSmoothSkygridLikelihood(String name, List<IntervalList> intervalList,
                                      Parameter logPopSizeParameter,
                                      Parameter gridPointParameter) {

        super(name);
        this.intervalList = intervalList;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;
        this.gridPoints = gridPointParameter.getParameterValues();
        this.numGridIntervals = gridPoints.length + 2;

        addVariable(logPopSizeParameter);
        addVariable(gridPointParameter);

        this.units = intervalList.get(0).getUnits();

        for (IntervalList intervals : intervalList) {
            if (intervals instanceof Model) {
                addModel((Model) intervals);
            }
            if (intervals.getUnits() != units) {
                throw new IllegalArgumentException("All intervalLists must have the same units.");
            }
        }

        if (!checkValidParameters(logPopSizeParameter, gridPointParameter)) {
            throw new IllegalArgumentException("Invalid initial parameters");
        }

        this.addKeyword("smooth0skygrid");
        if (intervalList.size() > 1) {
            this.addKeyword("multilocus");
        }
    }

    private List<TreeIntervals> debugIntervalList;
    public void setDebugIntervalList(List<TreeIntervals> intervals) {
        debugIntervalList = intervals;
    }

    public static boolean checkValidParameters(Parameter logPopSizes, Parameter gridPoints) {
        return (logPopSizes.getDimension() == gridPoints.getDimension() + 1) &&
                checkStrictlyIncreasing(gridPoints);
    }

    private static boolean checkStrictlyIncreasing(Parameter gridPoints) {
        double lastValue = gridPoints.getParameterValue(0);
        for (int index = 1; index < gridPoints.getDimension(); ++index) {
            double thisValue = gridPoints.getParameterValue(index);
            if (thisValue <= lastValue) {
                return false;
            }
            lastValue = thisValue;
        }
        return true;
    }

    @Override
    public Type getUnits() {
        return units;
    }

    @Override
    public void setUnits(Type units) {
        this.units = units;
        for (IntervalList intervals : intervalList) {
            intervals.setUnits(units);
        }
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }

    @Override
    protected double calculateLogLikelihood() {
        if (!intervalsKnown) {
            computeSufficientStatistics();
            intervalsKnown = true;
        }

        return 0;
    }

    class Event {
        double integratedIntensity;
        double logReciprocalPopSize;
        int lineages;
    }
    
    private void computeSufficientStatistics() {

        for (IntervalList intervals : intervalList) {

            // Find first grid-interval for tree
            double previousIntervalTimeOnTime = intervals.getStartTime();
            int currentGridIndex = 0;
            while (previousIntervalTimeOnTime > gridPoints[currentGridIndex]) {
                ++currentGridIndex;
            }

            // Interate over tree-intervals in time-order
            for (int j = 0; j < intervals.getIntervalCount() - 1; ++j) {
                double currentIntervalTimeOnTree = intervals.getIntervalTime(j);
                int currentIntervalLineageCount = intervals.getLineageCount(j);

                while(intervals.getIntervalTime(j + 1) > gridPoints[currentGridIndex]) {
                    // Do somethimg until end of internval
                }
                // Do something with last bit of time

                
            }

            // Handle last tree-interval here

        }

        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void storeState() {
        super.storeState();

        if (savedCoalescentCount == null) {
            savedCoalescentCount = new double[coalescentCount.length];
            savedCoalescentIntensity = new double[coalescentIntensity.length];
        }

        System.arraycopy(coalescentCount, 0, savedCoalescentCount, 0, coalescentCount.length);
        System.arraycopy(coalescentIntensity, 0, savedCoalescentIntensity, 0, coalescentIntensity.length);
    }

    @Override
    protected void restoreState() {
        super.restoreState();

        double[] tmp = coalescentCount;
        coalescentCount = savedCoalescentCount;
        savedCoalescentCount = tmp;

        tmp = coalescentIntensity;
        coalescentIntensity = savedCoalescentIntensity;
        savedCoalescentIntensity = tmp;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (variable == gridPointParameter) {
            throw new RuntimeException("Not yet implemented");
        }

        if (variable == logPopSizeParameter) {
            likelihoodKnown = false; // intervals are, however, still known
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);

        if (model instanceof IntervalList) {
            intervalsKnown = false;
            likelihoodKnown = false;
        }
    }

    @Override
    public int getNumberOfCoalescentEvents() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Differentiable skygrid coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(CommonCitations.GILL_2013_IMPROVING,
                new Citation(
                        new Author[] {
                                new Author("MA", "Suchard"),
                                new Author( "X", "Ji"),
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
    }

    public static double getIntensityInInterval(double time1, double time2,
                                                double startTime, double endTime,
                                                double startValue, double endValue,
                                                double beta) throws Exception {

        assert time1 >= startTime && time2 >= startTime;
        assert time1 <= endTime && time2 <= endTime;
        assert time1 != time2;

        if (beta == 1.0) {
            return getAnalyticIntensityForLinearModel(time1, time2, startTime, endTime, startValue, endValue);
        } else if (beta == 0.0) {
            return getAnalyticIntensityForConstantModel(time1, time2, startTime, endTime, startValue, endValue);
        } else if (beta == Double.POSITIVE_INFINITY) {
            return getAnalyticIntensityForShiftedSkygridModel(time1, time2, startTime, endTime, startValue, endValue);
        } else {
            return getNumericIntensityInInterval(time1, time2, startTime, endTime, startValue, endValue, beta);
        }
    }

    public static double getAnalyticIntensityForShiftedSkygridModel(double time1, double time2,
                                                                     double startTime, double endTime,
                                                                     double startValue, double endValue) {
        double midPoint = (startTime + endTime) * 0.5;
        if (time1 < midPoint) {
            if (time2 <= midPoint) {
                return (time2 - time1) / Math.exp(startValue);
            } else {
                return (midPoint - time1) / Math.exp(startValue) + (time2 - midPoint) / Math.exp(endValue);
            }
        } else {
            return (time2 - time1) / Math.exp(endValue);
        }
    }

    public static double getAnalyticIntensityForConstantModel(double time1, double time2,
                                                               double startTime, double endTime,
                                                               double startValue, double endValue) {
        return (time2 - time1) / Math.exp((startValue + endValue) * 0.5);
    }

    public static double[] getGradientWrtLogPopSizesInInterval(double time1, double time2,
                                                              double startTime, double endTime,
                                                              double startValue, double endValue,
                                                              double beta) throws Exception {
        double integratedFunctional = getNumericFunctionalInInterval(time1, time2,
                startTime, endTime, startValue, endValue, beta);

        double integratedIntensity = getNumericIntensityInInterval(time1, time2,
                startTime, endTime, startValue, endValue, beta);

        return new double[] { integratedFunctional - integratedIntensity, -integratedFunctional };
    }

    public static double[] getGradientWrtLogPopSizesInIntervalViaCentralDifference(double time1, double time2,
                                                                                   double startTime, double endTime,
                                                                                   double startValue, double endValue,
                                                                                   double beta) {

        double gradStartValue = NumericalDerivative.firstDerivative(new UnivariateFunction() {
            @Override
            public double evaluate(double x) {
                try {
                    return getIntensityInInterval(time1, time2, startTime, endTime, x, endValue, beta);
                } catch (Exception e) {
                    e.printStackTrace();
                }
                return 0;
            }

            @Override
            public double getLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
        }, startValue);

        double gradEndValue = NumericalDerivative.firstDerivative(new UnivariateFunction() {
            @Override
            public double evaluate(double x) {
                try {
                    return getIntensityInInterval(time1, time2, startTime, endTime, startValue, x, beta);
                } catch (Exception e) {
                    e.printStackTrace();
                }
                return 0;
            }

            @Override
            public double getLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
        }, endValue);

        return new double[] { gradStartValue, gradEndValue };
    }

    public static double getNumericIntensityInInterval(double time1, double time2,
                                                       double startTime, double endTime,
                                                       double startValue, double endValue,
                                                       double beta) throws Exception {
        UnivariateRealFunction f = v -> getReciprocalPopSizeInInterval(v, startTime, endTime,
                startValue, endValue, beta);
        return integrator.integrate(f, time1, time2);
    }

    public static double getNumericFunctionalInInterval(double time1, double time2,
                                                        double startTime, double endTime,
                                                        double startValue, double endValue,
                                                        double beta) throws Exception {
        UnivariateRealFunction f = v -> getReciprocalPopSizeInInterval(v, startTime, endTime,
                startValue, endValue, beta) / (1.0 + getScaledOdds(v, startTime, endTime, beta));
        return integrator.integrate(f, time1, time2);
    }

    public static double getNumericLogPopSizeIntensity(double time1, double time2,
                                                       double startTime, double endTime,
                                                       double beta) throws Exception {
        UnivariateRealFunction f = v -> getLogPopSizeIntensity(v, startTime, endTime, beta);
        return integrator.integrate(f, time1, time2);
    }

    public static double getAnalyticIntensityForLinearModel(double time1, double time2,
                                                            double startTime, double endTime,
                                                            double startValue, double endValue) {
        double slope = (endValue - startValue) / (endTime - startTime);
        time1 -= startTime;
        time2 -= startTime;
        double e1 = Math.exp(-slope * time1 - startValue);
        double e2 = Math.exp(-slope * time2 - startValue);
        return (e1 - e2) / slope;
    }

    public static double getReciprocalPopSizeInInterval(double time,
                                                        double startTime, double endTime,
                                                        double startValue, double endValue,
                                                        double beta) {
        return Math.exp(-getLogPopSizeInInterval(time, startTime, endTime, startValue, endValue, beta));
    }

    public static double getPopSizeInInterval(double time,
                                              double startTime, double endTime,
                                              double startValue, double endValue,
                                              double beta) {
        return Math.exp(getLogPopSizeInInterval(time, startTime, endTime, startValue, endValue, beta));
    }

    public static double getLogPopSizeIntensity(double time,
                                                double startTime, double endTime,
                                                double beta) {
        return Math.exp(-getLogPopSizeWeight(time, startTime, endTime, beta));
    }

    public static double getLogPopSizeWeight(double time,
                                             double startTime, double endTime,
                                             double beta) {
        assert time >= startTime;
        assert time <= endTime;

        if (time == startTime) { // Avoid divide-by-zero
            if (beta == 0.0) {
                return 2.0;
            } else {
                return Double.POSITIVE_INFINITY;
            }
        }

        return 1.0 + getScaledOdds(time, startTime, endTime, beta);
    }

    public static double getLogPopSizeInInterval(double time,
                                                 double startTime, double endTime,
                                                 double startValue, double endValue,
                                                 double beta) {
        assert time >= startTime;
        assert time <= endTime;

        if (time == startTime) { // Avoid divide-by-zero
            if (beta == 0.0) {
                return (startValue + endValue) * 0.5;
            } else {
                return startValue;
            }
        }

        double scaledOdds = getScaledOdds(time, startTime, endTime, beta);
        return startValue + (endValue - startValue) / (1.0 + scaledOdds);
    }

    private static double getScaledOdds(double time, double startTime, double endTime, double beta) {
        double timeOdds = (endTime - time) / (time - startTime);
        return Math.pow(timeOdds, beta);
    }
}
