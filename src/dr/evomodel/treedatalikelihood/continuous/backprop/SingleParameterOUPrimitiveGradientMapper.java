package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evolution.tree.NodeRef;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;

/**
 * Maps primitive OU gradients (dL/dS, dL/dΣ_stat, dL/dμ)
 * into the entries of a single BEAST Parameter.
 */
public interface SingleParameterOUPrimitiveGradientMapper {

    /** The BEAST Parameter this mapper writes into. */
    Parameter getParameter();

    /** Number of doubles this parameter contributes (usually parameter.getDimension()). */
    int getDimension();

    /**
     * Write this parameter's contribution into the global gradient vector.
     *
     * @param node       branch is parent(node)
     * @param dLdS       gradient w.r.t. S
     * @param dLdSigmaSt gradient w.r.t. Σ_stat
     * @param dLdMu      gradient w.r.t. μ
     * @param target     global gradient vector
     * @param offset     starting index in target for this parameter
     */
    void mapPrimitiveToParameter(NodeRef node,
                                 DenseMatrix64F dLdS,
                                 DenseMatrix64F dLdSigmaSt,
                                 DenseMatrix64F dLdMu,
                                 DenseMatrix64F dLdSigma,
                                 double[] target,
                                 int offset);
}
