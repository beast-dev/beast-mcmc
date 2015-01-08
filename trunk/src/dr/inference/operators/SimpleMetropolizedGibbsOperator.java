/**
 *
 */
package dr.inference.operators;

import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.prior.Prior;

import java.util.logging.Logger;

/**
 * @author Sebastian Hoehna
 */
public abstract class SimpleMetropolizedGibbsOperator extends SimpleOperator implements GeneralOperator {

    /**
     *
     */
    public SimpleMetropolizedGibbsOperator() {
        // Do nothing
    }

    public abstract double doOperation(Prior prior, Likelihood likelihood)
            throws OperatorFailedException;

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.GibbsOperator#getStepCount()
      */
    public abstract int getStepCount();

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.SimpleOperator#getOperatorName()
      */
    @Override
    public abstract String getOperatorName();

    public final double operate() throws OperatorFailedException {
        return operate(null, null);
    }

    public final double operate(Prior prior, Likelihood likelihood)
            throws OperatorFailedException {
        if (operateAllowed) {
            operateAllowed = false;
            return doOperation(prior, likelihood);
        } else
            throw new RuntimeException(
                    "Operate called twice without accept/reject in between!");
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMaximumAcceptanceLevel()
      */
    public final double getMaximumAcceptanceLevel() {
        return 1.0;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMaximumGoodAcceptanceLevel()
      */
    public final double getMaximumGoodAcceptanceLevel() {
        return 1.0;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMinimumAcceptanceLevel()
      */
    public final double getMinimumAcceptanceLevel() {
        return 0.005;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMinimumGoodAcceptanceLevel()
      */
    public final double getMinimumGoodAcceptanceLevel() {
        return 0.01;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
      */
    public final String getPerformanceSuggestion() {
        return "";
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getTargetAcceptanceProbability()
      */
    public final double getTargetAcceptanceProbability() {
        return 1.0;
    }

    protected double evaluate(Likelihood likelihood, Prior prior, double pathParameter) {

        double logPosterior = 0.0;

        if (prior != null) {
            final double logPrior = prior.getLogPrior(likelihood.getModel());

            if (logPrior == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }

            logPosterior += logPrior;
        }

        final double logLikelihood = likelihood.getLogLikelihood()*pathParameter;

        if (Double.isNaN(logLikelihood)) {
            return Double.NEGATIVE_INFINITY;
        }
        // System.err.println("** " + logPosterior + " + " + logLikelihood + " =
        // " + (logPosterior + logLikelihood));
        logPosterior += logLikelihood;

        return logPosterior;
    }

    protected void restore(Prior prior, Likelihood likelihood,
                           Model currentModel, MCMCOperator mcmcOperator, double oldScore) {
        currentModel.restoreModelState();

        // This is a test that the state is correctly restored. The restored
        // state is fully evaluated and the likelihood compared with that before
        // the operation was made.
        likelihood.makeDirty();
        final double testScore = evaluate(likelihood, prior, 1.0);

        if (Math.abs(testScore - oldScore) > 1e-6) {
            Logger.getLogger("error").severe(
                    "State was not correctly restored after reject step.\n"
                            + "Likelihood before: " + oldScore
                            + " Likelihood after: " + testScore + "\n"
                            + "Operator: " + mcmcOperator + " "
                            + mcmcOperator.getOperatorName());
        }
    }

}
