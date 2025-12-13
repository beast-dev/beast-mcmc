package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * Backpropagation-based diffusion gradient adapter for HMC.
 *
 * This class:
 *  - wraps a BranchSpecificGradient that internally uses ContinuousTraitBackpropGradient
 *  - exposes a single (possibly Compound) Parameter
 *  - returns the gradient slice corresponding to this parameter (using offset + dim)
 */
public final class BackPropParameterDiffusionGradient extends AbstractDiffusionGradient {

    private final BranchSpecificGradient branchSpecificGradient;

    /** The parameter this gradient provider is responsible for (may be CompoundParameter). */
    private final Parameter parameter;

    /** For reporting; usually the same as 'parameter'. */
    private final Parameter rawParameter;

    /** Dimension of 'parameter'. */
    private final int dim;

    private BackPropParameterDiffusionGradient(BranchSpecificGradient branchSpecificGradient,
                                               Likelihood likelihood,
                                               Parameter parameter,
                                               Parameter rawParameter,
                                               double upperBound,
                                               double lowerBound) {
        super(likelihood, upperBound, lowerBound);  // TODO vectorize lowerBound and upperBound
        this.branchSpecificGradient = branchSpecificGradient;

        this.parameter = parameter;
        this.rawParameter = rawParameter;
        this.dim = parameter.getDimension(); // TODO generalise
    }

    /**
     * Factory: creates a backprop diffusion gradient for a given Parameter
     * (often a CompoundParameter aggregating several primitives).
     */
    public static BackPropParameterDiffusionGradient createBackpropGradient(
            BranchSpecificGradient branchSpecificGradient,
            Likelihood likelihood,
            Parameter parameter) {

        double upper = Double.POSITIVE_INFINITY;
        double lower = Double.NEGATIVE_INFINITY;

        return new BackPropParameterDiffusionGradient(
                branchSpecificGradient,
                likelihood,
                parameter,
                parameter,
                upper,
                lower
        );
    }

    // ---------------------------------------------------------------------
    // GradientWrtParameterProvider / AbstractDiffusionGradient interface
    // ---------------------------------------------------------------------

    @Override
    public Parameter getParameter() {
        // HMC (and numeric checker) will read and update this parameter
        return parameter;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public Parameter getRawParameter() {
        return rawParameter;
    }


    /**
     * Public analytic gradient: returns this parameter's slice only.
     */
    @Override
    public double[] getGradientLogDensity() {
        // Full gradient for ALL parameters controlled by BranchSpecificGradient
        double[] full = branchSpecificGradient.getGradientLogDensity();
        return getGradientLogDensity(full);
    }

    /**
     * Slice this parameter's segment out of a full gradient vector using offset + dim.
     */
    @Override
    public double[] getGradientLogDensity(double[] gradient) {
        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = gradient[offset + i];
        }
        return result;
    }

    @Override
    public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter
    getDerivationParameter() {
        // For backprop we don't *really* use derivation parameters, but the interface
        // expects something. To minimise surprises, delegate if available.
        List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derList =
                branchSpecificGradient.getDerivationParameter();

        if (derList == null || derList.isEmpty()) {
            // Fallback dummy value; not used in the backprop computation itself.
            return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient
                    .DerivationParameter.WRT_GENERAL_SELECTION_STRENGTH;
        }
        return derList.get(0);
    }

    @Override
    public String getReport() {
        return "BackpropGradient." + rawParameter.getParameterName() + "\n" +
                super.getReport();
    }
}
