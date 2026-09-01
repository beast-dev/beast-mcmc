package dr.evomodel.epidemiology;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;

public class TwoPathogenCompatibleGrid extends Parameter.Abstract {

    private double[] gridStartTimes;
    private Bounds bounds = null;

    public TwoPathogenCompatibleGrid(Parameter cutOff,
                                     Parameter numGridPoints,
                                     Parameter mostRecentSamplingDateOne,
                                     Parameter mostRecentSamplingDateTwo) {

        this.gridStartTimes = computeGrid(cutOff.getParameterValue(0),
                numGridPoints.getParameterValue(0),
                mostRecentSamplingDateOne.getParameterValue(0),
                mostRecentSamplingDateTwo.getParameterValue(0));

    }

    protected double[] computeGrid(double cutOff, double numGridPoints, double mostRecentOne, double mostRecentTwo){
        double moreRecentTime = Math.max(mostRecentOne, mostRecentTwo);
        double lessRecentTime = Math.min(mostRecentOne, mostRecentTwo);
        double delta = cutOff/numGridPoints;
        double offset = moreRecentTime-lessRecentTime;
        int numSkippedIntervals = (int) Math.floor(offset/delta);
        int newNumGridPoints = (int) numGridPoints - numSkippedIntervals;
        double[] newGrid = new double[newNumGridPoints];

        newGrid[0] = 0;
        for (int i = 1; i < newNumGridPoints; i++) {
            newGrid[i] = (numSkippedIntervals + i)*delta - offset;
        }
        return newGrid;
    }

    protected void storeValues() {
    }

    protected void restoreValues() {
    }

    protected void acceptValues() {
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public int getDimension() {
        return gridStartTimes.length;
    }

    public double getParameterValue(int dim) {
        return gridStartTimes[dim];
    }

    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Not implemented");
    }

    public String getParameterName() {
        return "Grid points (starting times of intervals) in episodic birth-death sampling model for pathogen with narrower time scale in two-pathogen model";
    }

    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not implemented.");
    }
}
