package dr.inference.operators;

/**
 * Supplies the dependent-process edge evidence term from the supplement,
 *
 *   pre_b^T P_b(r) post_b,
 *
 * for CTMC dependent processes sharing the reward coordinate with the
 * modulating reward process.
 */
public interface RewardDependentCtmcEdgeEvidenceProvider {

    void prepare();

    double logEvidence(int branchNodeNumber, double rawReward);
}
