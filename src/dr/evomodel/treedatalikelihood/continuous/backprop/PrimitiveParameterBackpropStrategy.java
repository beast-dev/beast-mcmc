package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

public interface PrimitiveParameterBackpropStrategy {

    PrimitiveGradientSet backprop(
            ContinuousTraitBackpropGradient.OUBranchCache cache,
            MessageBackprop.Result leafGrads
    );

    final class PrimitiveGradientSet {
        public final DenseMatrix64F dLdS;
        public final DenseMatrix64F dLdSigma;
        public final DenseMatrix64F dLdSigmaStat;
        public final DenseMatrix64F dLdMu;

        public PrimitiveGradientSet(DenseMatrix64F dLdS, DenseMatrix64F dLdSigma, DenseMatrix64F dLdSigmaStat, DenseMatrix64F dLdMu) {
            this.dLdS = dLdS;
            this.dLdSigma = dLdSigma;
            this.dLdSigmaStat = dLdSigmaStat;
            this.dLdMu = dLdMu;
        }
    }
}

