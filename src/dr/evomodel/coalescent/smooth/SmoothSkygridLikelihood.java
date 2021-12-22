package dr.evomodel.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;

import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class SmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Citable {

    private final List<IntervalList> intervalList;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;
    private Type units;

    private final static RombergIntegrator integrator = new RombergIntegrator();
    private final static boolean USE_LINEAR_ANALYTIC_SOLUTION = true;

    public SmoothSkygridLikelihood(String name, List<IntervalList> intervalList,
                                   Parameter logPopSizeParameter,
                                   Parameter gridPointParameter) {

        super(name);
        this.intervalList = intervalList;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;

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

    @Override
    protected double calculateLogLikelihood() {
        return 0;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        if (variable == gridPointParameter) {
            throw new RuntimeException("Not yet implemented");
        }

        if (variable == logPopSizeParameter) {
            throw new RuntimeException("Not yet implemented");
        }

        throw new RuntimeException("Should not get here");
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model instanceof IntervalList) {
            intervalsKnown = false;
            likelihoodKnown = false;
        }

        throw new RuntimeException("Should not get here");
    }

    @Override
    public int getNumberOfCoalescentEvents() {
        return 0;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        return 0;
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        return null;
    }

    public static double getIntensityInInterval(double time1, double time2,
                                                double startTime, double endTime,
                                                double startValue, double endValue,
                                                double beta) throws Exception {

        if (USE_LINEAR_ANALYTIC_SOLUTION && beta == 1) {
            return getAnalyticIntensityForLinearModel(time1, time2, startTime, endTime, startValue, endValue);
        }

        UnivariateRealFunction f = v -> getReciprocalPopSizeInInterval(v, startTime, endTime,
                startValue, endValue, beta);

        return integrator.integrate(f, time1, time2);
    }

    private static double getAnalyticIntensityForLinearModel(double time1, double time2,
                                                             double startTime, double endTime,
                                                             double startValue, double endValue) {
        assert time1 >= startTime && time2 >= startTime;
        assert time1 <= endTime && time2 <= endTime;

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

    public static double getLogPopSizeInInterval(double time,
                                                 double startTime, double endTime,
                                                 double startValue, double endValue,
                                                 double beta) {
        assert time >= startTime;
        assert time <= endTime;

        if (time == startTime) { // Avoid divide-by-zero
            return startValue;
        }

        double scaledOdds = getScaledOdds(time, startTime, endTime, beta);
        return startValue + (endValue - startValue) / (1.0 + scaledOdds);
    }

    private static double getScaledOdds(double time, double startTime, double endTime, double beta) {
        double timeOdds = (endTime - time) / (time - startTime);
        return Math.pow(timeOdds, beta);
    }
}
