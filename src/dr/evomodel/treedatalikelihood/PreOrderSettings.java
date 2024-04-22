package dr.evomodel.treedatalikelihood;

/**
 * @author Marc A. Suchard
 */
public class PreOrderSettings {
    public boolean usePreOrder;
    public boolean branchRateDerivative;
    public boolean branchInfinitesimalDerivative;

    public boolean useAction;
    boolean useAmbiguities;


    public PreOrderSettings(boolean usePreOrder,
                            boolean branchRateDerivative,
                            boolean branchInfinitesimalDerivative,
                            boolean useAmbiguities,
                            boolean useAction) {
        this.usePreOrder = usePreOrder;
        this.branchRateDerivative = branchRateDerivative;
        this.branchInfinitesimalDerivative = branchInfinitesimalDerivative;
        this.useAmbiguities = useAmbiguities;
        this.useAction = useAction;
    }

    public static PreOrderSettings getDefault() {
        return new PreOrderSettings(false, false, false, false, false);
    }
}
