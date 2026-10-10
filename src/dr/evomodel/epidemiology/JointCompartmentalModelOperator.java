package dr.evomodel.epidemiology;

import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.JointOperator;

public class JointCompartmentalModelOperator extends JointOperator {

    private final CompartmentalModelSimulator simulator;

    public JointCompartmentalModelOperator (double weight, double targetAcceptanceProbability,
                                            CompartmentalModelSimulator simulator) {
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

    @Override
    public void addOperator(SimpleMCMCOperator operation) {
        operatorList.add(operation);
        if (operation instanceof AdaptableMCMCOperator &&
                ((AdaptableMCMCOperator) operation).getMode() != AdaptationMode.ADAPTATION_OFF) {
            operatorToOptimizeList.add(operatorList.size() - 1);
        }
    }

    @Override
    public AdaptationMode getMode() {
        boolean anyDefault = false;
        for (int i : operatorToOptimizeList) {
            AdaptationMode m = ((AdaptableMCMCOperator) operatorList.get(i)).getMode();
            if (m == AdaptationMode.ADAPTATION_ON) return AdaptationMode.ADAPTATION_ON;
            if (m == AdaptationMode.DEFAULT) anyDefault = true;
        }
        return anyDefault ? AdaptationMode.DEFAULT : AdaptationMode.ADAPTATION_OFF;
    }
}
