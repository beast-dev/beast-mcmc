package dr.evomodel.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class SmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Citable {

    private final IntervalList intervalList;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;

    public SmoothSkygridLikelihood(String name, IntervalList intervalList,
                                   Parameter logPopSizeParameter,
                                   Parameter gridPointParameter) {
        super(name, intervalList);
        this.intervalList = intervalList;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;

        addVariable(logPopSizeParameter);
        addVariable(gridPointParameter);
    }

    @Override
    public Type getUnits() {
        return intervalList.getUnits();
    }

    @Override
    public void setUnits(Type units) {
        intervalList.setUnits(units);
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
        super.handleVariableChangedEvent(variable, index, type);
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

    private static double getLogPopSizeInInterval(double time, double startTime, double endTime, double beta) {
        return 1.0 / getReciprocalLogPopSizeInInterval(time, startTime, endTime, beta);
    }

    private static double getLogPopSizeInInterval(double timeOdds, double beta) {
        return 1.0 / getReciprocalLogPopSizeInInterval(timeOdds, beta);
    }

    private static double getReciprocalLogPopSizeInInterval(double time, double startTime, double endTime, double beta) {
        assert time >= startTime;
        assert time < endTime;

        double timeOdds = (time - startTime) / (endTime - time);
        return getReciprocalLogPopSizeInInterval(timeOdds, beta);
    }

    private static double getReciprocalLogPopSizeInInterval(double timeOdds, double beta) {
        return 1.0 + Math.pow(timeOdds, beta);
    }
}
