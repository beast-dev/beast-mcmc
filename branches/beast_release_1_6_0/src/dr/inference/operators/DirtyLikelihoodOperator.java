package dr.inference.operators;

import dr.inference.model.Likelihood;
import dr.inferencexml.operators.DirtyLikelihoodOperatorParser;

/**
 * @author Marc Suchard
 */
public class DirtyLikelihoodOperator extends SimpleMCMCOperator implements GibbsOperator {

    public DirtyLikelihoodOperator(Likelihood likelihood, double weight) {
        this.likelihood = likelihood;
        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return DirtyLikelihoodOperatorParser.TOUCH_OPERATOR;
    }

    public double doOperation() throws OperatorFailedException {
        likelihood.makeDirty();
        return 0;
    }

    public int getStepCount() {
        return 1;
    }

    private Likelihood likelihood;
}
