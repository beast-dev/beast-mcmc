/*
 * NormalDistributionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.UnivariateFunction;
import dr.math.distributions.MarginalizedAlphaStableDistribution;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class MarginalizedAlphaStableDistributionModel extends AbstractModel
        implements ParametricDistributionModel, GradientProvider {

    private final Parameter scale;
    private final Parameter alpha;

    public MarginalizedAlphaStableDistributionModel(String name,
                                                    Parameter scale,
                                                    Parameter alpha) {
        super(name);
        this.scale = scale;
        this.alpha = alpha;

        addVariable(scale);
        addVariable(alpha);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    @Override
    public double logPdf(double[] x) {
        double logPdf = 0.0;
        for (double value : x) {
            logPdf += logPdf(value);
        }
        return logPdf;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    @Override
    public double logPdf(double x) {
        return MarginalizedAlphaStableDistribution.logPdf(x,
                scale.getParameterValue(0), alpha.getParameterValue(0));
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object obj) {
        double[] x;
        if (obj instanceof double[]) {
            x = (double[]) obj;
        } else {
            x = new double[1];
            x[0] = (Double) obj;
        }

        double[] result = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            result[i] = MarginalizedAlphaStableDistribution.gradLogPdf(x[i],
                    scale.getParameterValue(0), alpha.getParameterValue(0));
        }
        return result;
    }
    
    // TODO Notes for Gibbs sampling alpha:
//
//    if self.prior_type['global_shrinkage'] == 'jeffreys':
//
//        # Conjugate update for phi = 1 / gshrink ** reg_exponent
//        shape = beta_with_shrinkage.size / reg_exponent
//        scale = 1 / np.sum(np.abs(beta_with_shrinkage) ** reg_exponent)
//        phi = self.rg.np_random.gamma(shape, scale=scale)
//        gshrink = 1 / phi ** (1 / reg_exponent)
//
//    elif self.prior_type['global_shrinkage'] == 'half-cauchy':
//
//        gshrink = self.slice_sample_global_shrinkage(
//            gshrink, beta_with_shrinkage, self.prior_param['global_shrinkage']['scale'], reg_exponent
//        )
//



//    To update the global scale parameter τ , we work directly with the exponential-power density, marginalizing out the latent variables {ωj,uj}. This is a crucial source of efficiency in the bridge MCMC, and leads to the favorable mixing evident in Figure 1. From (1), observe that the posterior for ν ≡ τ−α, given β, is conditionally independent of y, and takes the form
//    p
//    p(ν | β) ∝ νp/α exp(−ν 􏰑|βj|α) p(ν).
//    j=1
//    Therefore if ν has a Gamma(c, d) prior, its conditional posterior will also be a gamma dis- tribution, with hyperparameters c⋆ = c+p/α and d⋆ = d+􏰏pj=1 |βj|α. To sample τ, simply draw ν from this gamma distribution, and use the transformation τ = ν−1/α.

//    def slice_sample_global_shrinkage(
//            self, gshrink, beta_with_shrinkage, global_scale, reg_exponent):
//        """ Slice sample phi = 1 / gshrink ** reg_exponent. """
//
//        n_update = 10 # Slice sample for multiple iterations to ensure good mixing.
//
//        # Initialize a gamma distribution object.
//        shape = (beta_with_shrinkage.size + 1) / reg_exponent
//        scale = 1 / np.sum(np.abs(beta_with_shrinkage) ** reg_exponent)
//        gamma_rv = sp.stats.gamma(shape, scale=scale)
//
//        phi = 1 / gshrink
//        for i in range(n_update):
//            u = self.rg.np_random.uniform() \
//                / (1 + (global_scale * phi ** (1 / reg_exponent)) ** 2)
//            upper = (np.sqrt(1 / u - 1) / global_scale) ** reg_exponent
//                # Invert the half-Cauchy density.
//            phi = gamma_rv.ppf(gamma_rv.cdf(upper) * self.rg.np_random.uniform())
//            if np.isnan(phi):
//                # Inverse CDF method can fail if the current conditional
//                # distribution is drastically different from the previous one.
//                # In this case, ignore the slicing variable and just sample from
//                # a Gamma.
//                phi = gamma_rv.rvs()
//        gshrink = 1 / phi ** (1 / reg_exponent)
//
//        return gshrink


//       def update_local_shrinkage(self, gshrink, beta_with_shrinkage, reg_exponent):
//
//        lshrink_sq = 1 / np.array([
//            2 * self.rg.tilted_stable(reg_exponent / 2, (beta_j / gshrink) ** 2)
//            for beta_j in beta_with_shrinkage
//        ])
//        lshrink = np.sqrt(lshrink_sq)
//
//        # TODO: Pick the lower and upper bound more carefully.
//        if np.any(lshrink == 0):
//            warn_message_only(
//                "Local shrinkage parameter under-flowed. Replacing with a small number.")
//            lshrink[lshrink == 0] = 10e-16
//        elif np.any(np.isinf(lshrink)):
//            warn_message_only(
//                "Local shrinkage parameter under-flowed. Replacing with a large number.")
//            lshrink[np.isinf(lshrink)] = 2.0 / gshrink
//
//        return lshrink
}
