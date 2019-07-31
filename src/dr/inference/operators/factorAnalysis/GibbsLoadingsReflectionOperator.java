package dr.inference.operators.factorAnalysis;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;

public class GibbsLoadingsReflectionOperator extends LoadingsReflectionOperator implements GibbsOperator {

    GibbsLoadingsReflectionOperator(SimpleMCMCOperator operator, MatrixParameterInterface loadings) {
        super(operator, loadings);
    }
}