package dr.evomodel.epidemiology;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;

public class TwoPathogenCompatibleNumGridPoints extends Parameter.Abstract {

    private double compatibleNumGridPoints;
    private Bounds bounds = null;

    public TwoPathogenCompatibleNumGridPoints(Parameter cutOff,
                                     Parameter numGridPoints,
                                     Parameter mostRecentSamplingDateOne,
                                     Parameter mostRecentSamplingDateTwo) {

        this.compatibleNumGridPoints = computeNumGridPoints(cutOff.getParameterValue(0),
                numGridPoints.getParameterValue(0),
                mostRecentSamplingDateOne.getParameterValue(0),
                mostRecentSamplingDateTwo.getParameterValue(0));

    }

    protected double computeNumGridPoints(double cutOff, double numGridPoints, double mostRecentOne, double mostRecentTwo){
        double moreRecentTime = Math.max(mostRecentOne, mostRecentTwo);
        double lessRecentTime = Math.min(mostRecentOne, mostRecentTwo);
        double delta = cutOff/numGridPoints;
        double offset = moreRecentTime-lessRecentTime;
        int numSkippedIntervals = (int) Math.floor(offset/delta);
        return numGridPoints - numSkippedIntervals;
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
        return 1;
    }

    public double getParameterValue(int dim) {
        return compatibleNumGridPoints;
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
        return "Number of grid points in episodic birth-death sampling model for pathogen with narrower time scale in two-pathogen model";
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
