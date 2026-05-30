package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

/**
 * System-property backed switches for dense fallback paths in block-diagonal
 * adjoint code.
 */
final class BlockDiagonalDenseFallbackPolicy {

    private static final String NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY =
            "beast.experimental.nativeForceDenseAdjointExp";
    private static final String TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY =
            "beast.experimental.transposeNativeFrechetInput";

    boolean forceDenseAdjointExp() {
        return Boolean.getBoolean(NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY);
    }

    boolean transposeNativeFrechetInput() {
        return Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY);
    }
}
