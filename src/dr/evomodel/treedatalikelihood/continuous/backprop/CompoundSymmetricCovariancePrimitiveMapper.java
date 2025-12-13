package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evolution.tree.NodeRef;
import dr.inference.model.CompoundParameter;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Primitive gradient mapper for a covariance matrix Σ that is
 * parameterized as a compound-symmetric correlation matrix:
 *
 *   Σ_ii = d_i
 *   Σ_ij = ρ * sqrt(d_i d_j),  i != j
 *
 * where d = diagonal parameter (length k), ρ = single off-diagonal parameter.
 *
 * This mapper consumes dL/dΣ from backprop and produces gradients
 * w.r.t. [d_0, ..., d_{k-1}, ρ].
 */
public final class CompoundSymmetricCovariancePrimitiveMapper
        implements SingleParameterOUPrimitiveGradientMapper {

    /** Underlying covariance/correlation matrix Σ (not the inverse). */
    private final CompoundSymmetricMatrix varianceMatrix;

    /** Diagonal variance parameter (length k). */
    private final Parameter diagParam;

    /** Single off-diagonal/correlation parameter ρ. */
    private final Parameter offDiagParam;

    /** Parameter exposed to HMC/backprop (diag and offDiag together). */
    private final CompoundParameter parameter;

    private final int dimTrait;   // k
    private final int dimParam;   // k + 1

    public CompoundSymmetricCovariancePrimitiveMapper(MatrixParameterInterface mi
//                                                      Parameter diagParam,
//                                                      Parameter offDiagParam
    ) {
        this.varianceMatrix = (CompoundSymmetricMatrix) mi;
        this.diagParam = varianceMatrix.getDiagonalParameter();//diagParam;
        this.offDiagParam = varianceMatrix.getUniqueParameter(1);//offDiagParam;

        this.dimTrait = diagParam.getDimension();
        this.dimParam = diagParam.getDimension() + offDiagParam.getDimension();

        // Expose [diag, offDiag] as a single parameter block
        this.parameter = new CompoundParameter(null);
        this.parameter.addParameter(diagParam);
        this.parameter.addParameter(offDiagParam);
    }

    @Override
    public Parameter getParameter() {
        // HMC sees [d_0,...,d_{k-1}, ρ] via this compound parameter
        return parameter;
    }

    @Override
    public int getDimension() {
        return dimParam;
    }

    @Override
    public void mapPrimitiveToParameter(NodeRef node,
                                        DenseMatrix64F dLdS,
                                        DenseMatrix64F dLdSigmaStat,
                                        DenseMatrix64F dLdMu,
                                        DenseMatrix64F dLdSigma, // we use THIS one
                                        double[] target,
                                        int offset) {

        // Safety: ensure dLdSigma is the right shape
        if (dLdSigma == null) {
            // No covariance gradient requested; fill zeros
            for (int i = 0; i < dimParam; i++) {
                target[offset + i] = 0.0;
            }
            return;
        }

        // Current diagonal and off-diagonal (ρ) values
        double[] d = diagParam.getParameterValues();
        double rho = offDiagParam.getParameterValue(0);

        // 1. Gradients w.r.t. diag entries d_i
        for (int i = 0; i < dimTrait; i++) {
            double gradDi = dLdSigma.get(i, i);  // contribution from Σ_ii = d_i

            double di = d[i];

            // Contributions from off-diagonal entries Σ_ij = ρ * sqrt(d_i d_j)
            for (int j = 0; j < dimTrait; j++) {
                if (j == i) continue;

                double dj = d[j];
                double g_ij = dLdSigma.get(i, j);
                double g_ji = dLdSigma.get(j, i);

                // ∂Σ_ij/∂d_i = ρ * (1/2) * sqrt(d_j / d_i)
                double dSigma_ij_dDi = 0.5 * rho * Math.sqrt(dj / di);

                // total contribution: g_ij * ∂Σ_ij/∂d_i + g_ji * ∂Σ_ji/∂d_i
                gradDi += (g_ij + g_ji) * dSigma_ij_dDi;
            }

            target[offset + i] = gradDi;
        }

        // 2. Gradient w.r.t. ρ (offDiag parameter)
        double gradRho = 0.0;
        for (int i = 0; i < dimTrait; i++) {
            double di = d[i];
            for (int j = i + 1; j < dimTrait; j++) {
                double dj = d[j];
                double g_ij = dLdSigma.get(i, j);
                double g_ji = dLdSigma.get(j, i);

                double sqrtDidj = Math.sqrt(di * dj);

                // ∂Σ_ij/∂ρ = sqrt(d_i d_j)
                // contributions from (i,j) and (j,i)
                gradRho += (g_ij + g_ji) * sqrtDidj;
            }
        }

        target[offset + dimTrait] = gradRho;
    }
}
