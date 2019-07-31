package dr.inference.operators.factorAnalysis;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.GeneralOperator;
import dr.inference.operators.SimpleMCMCOperator;

public class GeneralLoadingsReflectionOperator extends LoadingsReflectionOperator implements GeneralOperator {

    GeneralLoadingsReflectionOperator(SimpleMCMCOperator operator, MatrixParameterInterface loadings) {
        super(operator, loadings);
    }
}