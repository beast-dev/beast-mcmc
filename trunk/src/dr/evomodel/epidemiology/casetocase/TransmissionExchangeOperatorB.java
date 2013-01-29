package dr.evomodel.epidemiology.casetocase;

import dr.evomodel.operators.AbstractTreeOperator;
import dr.inference.operators.OperatorFailedException;

/**
 * Implements branch exchange operations that also exchange entire subtrees of the transmission tree. As this already
 * severely restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange.
 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorB extends AbstractTreeOperator {

    private final CaseToCaseTransmissionLikelihood c2cLikelihood;

    public TransmissionExchangeOperatorB(CaseToCaseTransmissionLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public void wide() throws OperatorFailedException {

    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
