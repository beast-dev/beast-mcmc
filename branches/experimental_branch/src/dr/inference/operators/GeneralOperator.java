package dr.inference.operators;

import dr.inference.model.Likelihood;
import dr.inference.prior.Prior;

/**
 * An operator that requires both the prior and likelihood to calculate hasting ratio.
 *
 * @author Alexei Drummond
 */
public interface GeneralOperator {

    double operate(Prior prior, Likelihood likelihood) throws OperatorFailedException;
}
