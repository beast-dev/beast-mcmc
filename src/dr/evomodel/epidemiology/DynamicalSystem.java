package dr.evomodel.epidemiology;

import java.util.ArrayList;
import java.util.List;

public class DynamicalSystem {

    private List<DynamicalVariable> variables = new ArrayList<DynamicalVariable>();
    private List<DynamicalForce> forces = new ArrayList<DynamicalForce>();
    private double currentTime = 0.0;
    private double timeStep = 0.0;

    public static void main(String[] args) {

        DynamicalSystem syst = new DynamicalSystem(0, 0.1);

        double time = syst.getTime();
        double beta = 0.25;
        double gamma = 0.2;

        DynamicalVariable susceptibles = new DynamicalVariable("susceptibles", time, 999.0);
        syst.addVariable(susceptibles);

        DynamicalVariable infecteds = new DynamicalVariable("infecteds", time, 1.0);
        syst.addVariable(infecteds);

        DynamicalVariable recovereds = new DynamicalVariable("recovereds", time, 0.0);
        syst.addVariable(recovereds);

        DynamicalVariable total = new DynamicalVariable("total", time, 1000.0);
        syst.addVariable(total);

        DynamicalForce f1 = new DynamicalForce("f1", -1*beta, susceptibles);
        f1.addMultiplier(infecteds);
        f1.addMultiplier(susceptibles);
        f1.addDivisor(total);
        syst.addForce(f1);

        DynamicalForce f2 = new DynamicalForce("f2", beta, infecteds);
        f2.addMultiplier(infecteds);
        f2.addMultiplier(susceptibles);
        f2.addDivisor(total);
        syst.addForce(f2);

        DynamicalForce f3 = new DynamicalForce("f3", -1*gamma, infecteds);
        f3.addMultiplier(infecteds);
        syst.addForce(f3);

        DynamicalForce f4 = new DynamicalForce("f4", gamma, recovereds);
        f4.addMultiplier(infecteds);
        syst.addForce(f4);

        for (int i=0; i<10; i++) {
            syst.step();
        }
        syst.print(0,1,0.1);

    }

    public DynamicalSystem(double t, double dt) {
        currentTime = t;
        timeStep = dt;
    }

    public double getTime() {
        return currentTime;
    }

    public void addVariable(DynamicalVariable var) {
        variables.add(var);
    }

    public void addForce(DynamicalForce frc) {
        forces.add(frc);
    }

    public void step() {
        for (DynamicalForce frc : forces) {
            frc.modCurrentValue(currentTime, timeStep);
        }
        for (DynamicalVariable var : variables) {
            var.modCurrentTime(timeStep);
            var.pushCurrentState();
        }
        currentTime += timeStep;
    }

    public void print(double start, double finish, double step) {

        System.out.print("time");
        for (DynamicalVariable var : variables) {
            System.out.print("\t" + var.getName());
        }
        System.out.println();

        for (double t=start; t<=finish; t += step) {
            System.out.printf("%.3f", t);
            for (DynamicalVariable var : variables) {
                System.out.printf("\t%.3f", var.getValue(t));
            }
            System.out.println();
        }

    }

}
