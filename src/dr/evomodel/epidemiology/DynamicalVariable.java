package dr.evomodel.epidemiology;

import java.util.ArrayList;
import java.util.List;
import java.util.Collections;

public class DynamicalVariable {

    private String name = "variable";
    private List<Double> times = new ArrayList<Double>();
    private List<Double> values = new ArrayList<Double>();
    private double currentTime = 0.0;
    private double currentValue = 0.0;

    // initialize variable
    public DynamicalVariable(String n, double t0, double v0) {
        name = n;
        currentTime = t0;
        currentValue = v0;
        add(t0, v0);
    }

    // add a time / value pair
    public void add(double t, double v) {
        Double tD = new Double(t);
        Double vD = new Double(v);
        times.add(tD);
        values.add(vD);
    }

    public void pushCurrentState() {
        add(currentTime,currentValue);
    }

    public void modCurrentTime(double dt) {
        currentTime += dt;
    }

    public void modCurrentValue(double dv) {
        currentValue += dv;
    }

    // reset arrays
    public void reset(double t0, double v0) {
        times.clear();
        values.clear();
        add(t0, v0);
    }

    public int size() {
        return times.size();
    }

    public String getName() {
        return name;
    }

    public double getTimeAtIndex(int i) {
        return times.get(i).doubleValue();
    }

    public double getValueAtIndex(int i) {
        return values.get(i).doubleValue();
    }

    // pull value closest to given point in time
    // uses linear interpolation when the exact time point isn't present
    public double getValue(double t) {

        double rval = 0.0;
        Double tD = new Double(t);
        int index = Collections.binarySearch(times, tD);

        // found exact match
        if (index >= 0) {
            Double vD = values.get(index);
            rval = vD.doubleValue();
        }

        // need to interpolate between two adjacent values
        else {
            int a = Math.abs(index) - 2;
            int b = Math.abs(index) - 1;
            double t0 = getTimeAtIndex(a);
            double v0 = getValueAtIndex(a);
            double t1 = getTimeAtIndex(b);
            double v1 = getValueAtIndex(b);
            rval = v0 + (t-t0) * ( (v1-v0) / (t1-t0) );
        }

        return rval;
    }

    // return integral between two points in time
    // approximation via trapezoidal rule
    public double getIntegral(double start, double finish) {

        // first index after start and last index before finish
        int a = 0;
        int b = 0;
        double sum = 0.0;

        // if not exact match, need first index after than start
        int index = Collections.binarySearch(times, new Double(start));
        if (index >= 0)
            a = index;
        else
            a = Math.abs(index) - 1;

         // if not exact match, need last index before than finish
        index = Collections.binarySearch(times, new Double(finish));
        if (index >= 0)
            b = index;
        else
            b = Math.abs(index) - 2;

        // from start to a
        sum += 0.5 * (getTimeAtIndex(a) - start) * (getValueAtIndex(a) + getValue(start));

        // between a and b
        for (index = a; index < b; index++) {
            sum += 0.5 * (getTimeAtIndex(index+1) - getTimeAtIndex(index)) * (getValueAtIndex(index+1) + getValueAtIndex(index));
        }

        // b to finish
        sum += 0.5 * (finish - getTimeAtIndex(b)) * (getValue(finish) + getValueAtIndex(b));

        return sum;

    }

    public void print() {
        System.out.println(name);
        for (int i=0; i<size(); i++) {
            System.out.println(times.get(i) + "\t" + values.get(i));
        }
    }

}
