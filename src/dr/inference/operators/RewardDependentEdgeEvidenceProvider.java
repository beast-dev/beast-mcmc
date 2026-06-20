package dr.inference.operators;

/**
 * Supplies a dependent-process evidence term for a candidate raw reward on one
 * branch. Implementations may evaluate the term locally or by exact likelihood
 * re-evaluation, but constants common to all candidates for the same branch are
 * allowed because the Gibbs weights are normalized over candidates.
 */
public interface RewardDependentEdgeEvidenceProvider {

    void prepare();

    double logEvidence(int branchNodeNumber, double rawReward);
}
