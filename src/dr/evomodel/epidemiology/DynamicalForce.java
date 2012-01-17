package dr.evomodel.epidemiology;

import java.util.ArrayList;
import java.util.List;

public class DynamicalForce {

    private String name = "force";
    private double coefficient = 0.0;
    private List<DynamicalVariable> multipliers = new ArrayList<DynamicalVariable>();
    private List<DynamicalVariable> divisors = new ArrayList<DynamicalVariable>();
    private DynamicalVariable targetVariable;

    public DynamicalForce(String n, double c, DynamicalVariable target) {
        name = n;
        coefficient = c;
        targetVariable = target;
    }

    public String getName() {
        return name;
    }

    public void addMultiplier(DynamicalVariable var) {
        multipliers.add(var);
    }

    public void addDivisor(DynamicalVariable var) {
        divisors.add(var);
    }

    public double getForce(double t) {
        double force = coefficient;
        for (DynamicalVariable var : multipliers) {
            force *= var.getValue(t);
        }
        for (DynamicalVariable var : divisors) {
            force /= var.getValue(t);
        }
        return force;
    }

    public void modCurrentValue(double t, double dt) {
        double dv = getForce(t) * dt;
        targetVariable.modCurrentValue(dv);

    }

}
