package dr.evomodel.epidemiology.casetocase;

import dr.evomodel.operators.AbstractTreeOperator;
import dr.inference.operators.OperatorFailedException;

/**
 * Implements the Wilson-Balding branch swapping move if the relocated phylogenetic subtree corresponds to a
 * transmission subtree.
 *
 * @author Matthew Hall
 */

public class TransmissionWilsonBaldingB extends AbstractTreeOperator{
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
