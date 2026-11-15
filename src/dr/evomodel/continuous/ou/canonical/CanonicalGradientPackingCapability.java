package dr.evomodel.continuous.ou.canonical;

/**
 * Capability for native canonical selection-gradient dimensions and final
 * packing into the requested parameter layout.
 */
public interface CanonicalGradientPackingCapability {

    int getSelectionGradientDimension();

    int getCompressedSelectionGradientDimension();

    int getNativeSelectionGradientScratchDimension();

    void finishNativeSelectionGradient(double[] compressedGradient,
                                       double[] nativeGradientScratch,
                                       double[] rotationGradientFlat,
                                       double[] gradientOut);
}
