/*
 * RewardDensityDomainException.java
 *
 * Signals that a reward-mixture branch model reached a mathematically
 * invalid input for the Sericola reward-density calculation (a reward
 * value that is NaN, or out of [0,1] by more than ordinary floating-point
 * noise). This is treated the same way BEAST's other HMC operators treat
 * numerical instability during a trajectory (see
 * HamiltonianMonteCarloOperator.NumericInstabilityException): the current
 * proposal is rejected, not the whole chain crashed. Unchecked so it can
 * propagate through the tree-likelihood call stack (TreeDataLikelihood,
 * DiscreteDataLikelihoodDelegate, RewardsAwarePartialsRepresentation, ...)
 * without requiring throws declarations on every intermediate method.
 */

package dr.inference.markovjumps;

public class RewardDensityDomainException extends RuntimeException {

    public RewardDensityDomainException(String message) {
        super(message);
    }
}
