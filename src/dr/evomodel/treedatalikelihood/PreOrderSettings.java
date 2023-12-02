package dr.evomodel.treedatalikelihood;

/**
 * @author Marc A. Suchard
 */
public class PreOrderSettings {
    boolean usePreOrder;
    boolean branchRateDerivative;
    boolean branchInfinitesimalDerivative;

    boolean useAmbiguities;

    public PreOrderSettings(boolean usePreOrder,
                            boolean branchRateDerivative,
                            boolean branchInfinitesimalDerivative,
                            boolean useAmbiguities) {
        this.usePreOrder = usePreOrder;
        this.branchRateDerivative = branchRateDerivative;
        this.branchInfinitesimalDerivative = branchInfinitesimalDerivative;
        this.useAmbiguities = useAmbiguities;
    }

    public static PreOrderSettings getDefault() {
        return new PreOrderSettings(false, false, false, false);
    }
}
