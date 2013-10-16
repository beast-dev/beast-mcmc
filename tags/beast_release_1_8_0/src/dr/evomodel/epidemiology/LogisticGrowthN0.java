/*
 * LogisticGrowthN0.java
 *
 * Daniel Wilson 4th October 2011
 *
 */

package dr.evomodel.epidemiology;

import dr.evolution.coalescent.*;

/**
 * This class models logistic growth.
 *
 * @author Daniel Wilson
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: LogisticGrowth.java,v 1.15 2008/03/21 20:25:56 rambaut Exp $
 */
public class LogisticGrowthN0 extends ExponentialGrowth {

    /**
     * Construct demographic model with default settings
     */
    public LogisticGrowthN0(Type units) {
        super(units);
    }

	public void setT50(double value) {
		t50 = value;
	}
	
	public double getT50() {
		return t50;
	}

    // Implementation of abstract methods

    /**
     * Gets the value of the demographic function N(t) at time t.
     *
     * @param t the time
     * @return the value of the demographic function N(t) at time t.
     */
    public double getDemographic(double t) {
        double N0 = getN0();
        double r = getGrowthRate();
        double T50 = getT50();
		
		return N0 * (1 + Math.exp(-r * T50)) / (1 + Math.exp(r * (t-T50)));
    }

    public double getLogDemographic(double t) {
		return Math.log(getDemographic(t));
	}

    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {
        double N0 = getN0();
        double r = getGrowthRate();
        double T50 = getT50();
		double exp_rT50 = Math.exp(-r*T50);
		
		return (t + exp_rT50 * (Math.exp(r * t) - 1)/r) / (N0 * (1 + exp_rT50));
    }

    /**
     * Returns the inverse function of getIntensity
     */
    public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public double getIntegral(double start, double finish) {
		return getIntensity(finish) - getIntensity(start);
    }

    public int getNumArguments() {
        return 3;
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0:
                return "N0";
            case 1:
                return "r";
            case 2:
                return "t50";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN0();
            case 1:
                return getGrowthRate();
            case 2:
                return getT50();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public void setArgument(int n, double value) {
        switch (n) {
            case 0:
                setN0(value);
                break;
            case 1:
                setGrowthRate(value);
                break;
            case 2:
                setT50(value);
                break;
            default:
                throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    public DemographicFunction getCopy() {
        LogisticGrowthN0 df = new LogisticGrowthN0(getUnits());
        df.setN0(getN0());
        df.setGrowthRate(getGrowthRate());
        df.setT50(getT50());
        return df;
    }

    //
    // private stuff
    //

    private double t50;
}
