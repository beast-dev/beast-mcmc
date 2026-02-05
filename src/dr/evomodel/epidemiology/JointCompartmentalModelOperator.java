package dr.evomodel.epidemiology;

import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.JointOperator;

public class JointCompartmentalModelOperator extends JointOperator {

    CompartmentalModel compartmentalModel;

    public JointCompartmentalModelOperator (double weight, double targetAcceptanceProbability,
                                            CompartmentalModel compartmentalModel) {
        super(weight, targetAcceptanceProbability);

        //operatorList = new ArrayList<SimpleMCMCOperator>();
        //operatorToOptimizeList = new ArrayList<Integer>();

        this.compartmentalModel = compartmentalModel;

        setWeight(weight);
    }

    public double doOperation() {

        compartmentalModel.simulateTrajectory();

        double logP = 0;

        for (SimpleMCMCOperator operation : operatorList) {

            logP += operation.doOperation();
        }

        return logP;
    }
}
