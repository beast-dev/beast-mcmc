package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evomodel.continuous.MultivariateElasticModel;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Declarative branch-to-OU-process specification for canonical OU likelihoods.
 *
 * <p>The first supported heterogeneous form keeps diffusion and stationary mean
 * shared, and varies only the selection model by branch. This mirrors the current
 * gradient contract while giving the tree likelihood a real branch-process map.
 */
public final class CanonicalOUBranchProcessSpecification {

    private final List<MultivariateElasticModel> elasticModels;
    private final Parameter branchProcessIndexParameter;

    public CanonicalOUBranchProcessSpecification(
            final List<MultivariateElasticModel> elasticModels,
            final Parameter branchProcessIndexParameter) {
        if (elasticModels == null || elasticModels.isEmpty()) {
            throw new IllegalArgumentException("At least one OU branch process is required.");
        }
        if (branchProcessIndexParameter == null) {
            throw new IllegalArgumentException("Branch-process index parameter is required.");
        }
        this.elasticModels =
                Collections.unmodifiableList(new ArrayList<MultivariateElasticModel>(elasticModels));
        this.branchProcessIndexParameter = branchProcessIndexParameter;
    }

    public int getProcessCount() {
        return elasticModels.size();
    }

    public MultivariateElasticModel getElasticModel(final int processIndex) {
        return elasticModels.get(processIndex);
    }

    public List<MultivariateElasticModel> getElasticModels() {
        return elasticModels;
    }

    public Parameter getBranchProcessIndexParameter() {
        return branchProcessIndexParameter;
    }
}
