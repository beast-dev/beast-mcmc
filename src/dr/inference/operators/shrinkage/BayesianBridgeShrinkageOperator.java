package dr.inference.operators.shrinkage;

import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;

/**
 * @author Marc A. Suchard
 */

public class BayesianBridgeShrinkageOperator extends SimpleMCMCOperator implements GibbsOperator {

    @Override
    public String getOperatorName() {
        return null;
    }

    @Override
    public double doOperation() {
        return 0;
    }
}
