package dr.evomodel.treedatalikelihood;

/**
 * @author Marc A. Suchard
 */
public class PreOrderSettings {
    boolean usePreOrder;
    boolean branchRateDerivative;
    boolean branchInfinitesimalDerivative;

    boolean useAction;

    public PreOrderSettings(boolean usePreOrder, boolean branchRateDerivative, boolean branchInfinitesimalDerivative,
                            boolean useAction) {
        this.usePreOrder = usePreOrder;
        this.branchRateDerivative = branchRateDerivative;
        this.branchInfinitesimalDerivative = branchInfinitesimalDerivative;
    }

    public static PreOrderSettings getDefault() {
        return new PreOrderSettings(false, false, false, false);
    }
}
