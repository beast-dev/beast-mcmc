package dr.evomodel.epidemiology.casetocase;

import dr.evolution.coalescent.DemographicFunction;
import dr.evomodel.epidemiology.LogisticGrowthN0;

/**
 * Shifted, four-parameter logistic growth function.
 *
 * r, growth rate parameter
 * t50, middle of the exponential phase (intended to be the start of the process).
 * n50, effective population size at t50
 * n0, carrying capacity
 */

public class LogisticGrowthN0N50 extends LogisticGrowthN0 {

    private double n50;

    public LogisticGrowthN0N50(Type units) {
        super(units);
    }

    public double getN50() {
        return n50;
    }

    public void setN50(double n50) {
        this.n50 = n50;
    }

    public double getDemographic(double t) {
        double r = getGrowthRate();
        double t50 = getT50();
        double n0 = getN0();

        double answer = 2*(n0 - n50)/(1+Math.exp(-r*(t50-t))) - (n0 - 2*n50);

        return answer;
    }


    public double getLogDemographic(double t) {
        double r = getGrowthRate();
        double t50 = getT50();
        double n0 = getN0();

        return Math.log(n0+(n0-2*n50)*(1+Math.exp(-r*(t50-t)))) - Math.log(1+Math.exp(-r*(t50-t)));
    }

    // warning: these are likely to give numerical errors if (N0-2N50)*exp(-r*T50) + N0 <= 0. The model class is
    // parametrised to avoid this ever happening

    public double getIntensity(double t) {

        double r = getGrowthRate();
        double t50 = getT50();
        double n0 = getN0();

        return (t/n0)
                +(2*n50/(r*n0*(n0-2*n50)))
                *Math.log((((n0-2*n50)*Math.exp(-r*(t50-t)))+n0)
                /(((n0-2*n50)*Math.exp(-r*t50))+n0));
    }

    public double getIntegral(double start, double finish) {

        double r = getGrowthRate();
        double t50 = getT50();
        double n0 = getN0();

        return ((finish-start)/n0)
                +(2*n50/(r*n0*(n0-2*n50)))
                *Math.log((((n0-2*n50)*Math.exp(-r*(t50-finish)))+n0)
                /(((n0-2*n50)*Math.exp(-r*(t50-start)))+n0));

    }


    public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public int getNumArguments() {
        return 4;
    }

    public String getArgumentName(int n) {
        switch (n) {
            case 0:
                return "N_T";
            case 1:
                return "r";
            case 2:
                return "t50";
            case 3:
                return "N_limit";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public double getArgument(int n) {
        switch (n) {
            case 0:
                return getN50();
            case 1:
                return getGrowthRate();
            case 2:
                return getT50();
            case 3:
                return getN0();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    public void setArgument(int n, double value) {
        switch (n) {
            case 0:
                setN50(value);
                break;
            case 1:
                setGrowthRate(value);
                break;
            case 2:
                setT50(value);
                break;
            case 3:
                setN0(value);
                break;
            default:
                throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        //The process started at N_T
        return getN0();
    }

    public DemographicFunction getCopy() {
        LogisticGrowthN0N50 df = new LogisticGrowthN0N50(getUnits());
        df.setN50(getN50());
        df.setGrowthRate(getGrowthRate());
        df.setT50(getT50());
        df.setN0(getN0());
        return df;
    }

}
