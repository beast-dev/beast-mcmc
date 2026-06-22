package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;

/**
 * Branch-local view of tree messages needed to build canonical transition adjoints.
 */
public final class BranchGradientMessageBundle {

    private int childNodeIndex;
    private boolean externalChild;
    private boolean fixedParent;
    private CanonicalGaussianState abovePreOrderAtParent;
    private CanonicalGaussianState belowPostOrderAtChild;
    private CanonicalTipObservationModel tipObservationModel;
    private double[] fixedRootValue;

    public void set(final int childNodeIndex,
                    final boolean externalChild,
                    final boolean fixedParent,
                    final CanonicalGaussianState abovePreOrderAtParent,
                    final CanonicalGaussianState belowPostOrderAtChild,
                    final CanonicalTipObservationModel tipObservationModel,
                    final double[] fixedRootValue) {
        this.childNodeIndex = childNodeIndex;
        this.externalChild = externalChild;
        this.fixedParent = fixedParent;
        this.abovePreOrderAtParent = abovePreOrderAtParent;
        this.belowPostOrderAtChild = belowPostOrderAtChild;
        this.tipObservationModel = tipObservationModel;
        this.fixedRootValue = fixedRootValue;
    }

    public int getChildNodeIndex() {
        return childNodeIndex;
    }

    public boolean isExternalChild() {
        return externalChild;
    }

    public boolean isFixedParent() {
        return fixedParent;
    }

    public CanonicalGaussianState getAbovePreOrderAtParent() {
        return abovePreOrderAtParent;
    }

    public CanonicalGaussianState getBelowPostOrderAtChild() {
        if (belowPostOrderAtChild == null) {
            throw new IllegalStateException(
                    "Branch " + childNodeIndex + " does not have a finite below-message state.");
        }
        return belowPostOrderAtChild;
    }

    public CanonicalTipObservationModel getTipObservationModel() {
        if (tipObservationModel == null) {
            throw new IllegalStateException(
                    "Branch " + childNodeIndex + " does not have a tip observation model.");
        }
        return tipObservationModel;
    }

    public double[] getFixedRootValue() {
        if (!fixedParent || fixedRootValue == null) {
            throw new IllegalStateException(
                    "Branch " + childNodeIndex + " does not have a fixed parent value.");
        }
        return fixedRootValue;
    }
}
