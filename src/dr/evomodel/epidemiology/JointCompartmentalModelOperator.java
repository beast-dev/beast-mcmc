package dr.evomodel.epidemiology;

import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.JointOperator;

public class JointCompartmentalModelOperator extends JointOperator {

    //CompartmentalModel compartmentalModel;
    StochasticSimulator simulator;

    public JointCompartmentalModelOperator (double weight, double targetAcceptanceProbability,
                                            StochasticSimulator simulator) {
        super(weight, targetAcceptanceProbability);

        //operatorList = new ArrayList<SimpleMCMCOperator>();
        //operatorToOptimizeList = new ArrayList<Integer>();

        this.simulator = simulator;

        setWeight(weight);
    }

    public double doOperation() {

        simulator.simulateTrajectory();

        double logP = 0;

        for (SimpleMCMCOperator operation : operatorList) {

            logP += operation.doOperation();
        }

        return logP;
    }
}
