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

    // pull value closest to given point in time
    public double getValue(double t) {
        Double tD = new Double(t);
        int index = Collections.binarySearch(times, tD);
        Double vD = values.get(index);
        return vD.doubleValue();
    }

    public void print() {
        System.out.println(name);
        for (int i=0; i<size(); i++) {
            System.out.println(times.get(i) + "\t" + values.get(i));
        }
    }

}
