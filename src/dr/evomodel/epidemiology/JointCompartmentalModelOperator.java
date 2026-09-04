package dr.evomodel.epidemiology;

import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.JointOperator;

public class JointCompartmentalModelOperator extends JointOperator {

    private final StochasticSimulator simulator;

    public JointCompartmentalModelOperator (double weight, double targetAcceptanceProbability,
                                            StochasticSimulator simulator) {
        super(weight, targetAcceptanceProbability);
        this.simulator = simulator;
        setWeight(weight);
    }

    public double doOperation() {
        double logP = 0;
        for (SimpleMCMCOperator operation : operatorList) {

            logP += operation.doOperation();
        }
        simulator.simulateTrajectory();
        return logP;
    }
}
