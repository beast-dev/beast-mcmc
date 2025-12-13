package dr.evomodel.continuous;

import org.ejml.data.DenseMatrix64F;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * Immutable value object representing a basis decomposition.
 * For block-diagonal: A = R D R^{-1}, with D provided in compressed form.
 */
/*
    * @author Filippo Monti
 */
public final class BasisRepresentation {

    private final int dim;
    private final int nParameters;

    // For block-diagonal: valuesD is compressed [diag | upper | lower], length = 3*dim-2
    // For non-block strategies: you can still store "eigenvalues" in valuesD (length dim), if you want,
    // but then blockStructure must be null and consumers must interpret accordingly.
    private final double[] valuesD;

    // Always row-major dim*dim
    private final double[] valuesR;
    private final double[] valuesRinv; // can be null for strategies that do not require it

    private final BlockStructure blockStructure; // null if not block-diagonal

    public BasisRepresentation(int dim,
                               int nParameters,
                               double[] valuesD,
                               double[] valuesR,
                               double[] valuesRinv) {
        this(dim, nParameters, valuesD, valuesR, valuesRinv, null);
    }

    public BasisRepresentation(int dim,
                               int nParameters,
                               double[] valuesD,
                               double[] valuesR,
                               double[] valuesRinv,
                               BlockStructure blockStructure) {
        this.dim = dim;
        this.nParameters = nParameters;
        // TODO CHECK THIS CHANGE
        this.valuesD = valuesD;
        this.valuesR = valuesR;
        this.valuesRinv = valuesRinv;
//        this.valuesD = valuesD.clone();
//        this.valuesR = (valuesR == null ? null : valuesR.clone());
//        this.valuesRinv = (valuesRinv == null ? null : valuesRinv.clone());
        this.blockStructure = blockStructure;
    }

    // ---- Legacy wrappers (avoid if you can; these wrap the underlying arrays) ----

    public DenseMatrix64F getR() {
        return wrap(valuesR, 0, dim, dim);
    }

    public DenseMatrix64F getRinv() {
        if (valuesRinv == null) {
            return null;
        }
        return wrap(valuesRinv, 0, dim, dim);
    }

    // ---- Accessors ----

    /** Legacy name: for block-diagonal this is NOT eigenvalues; it is valuesD as stored. */
    public double[] getEigenvalues() {
        return valuesD;
    }

    public double[] getEigenvectors() {
        return valuesR;
    }

    public int getNParameters() {
        return nParameters;
    }

    /**
     * Returns:
     * - if block-diagonal: concatenation [R | Rinv], both dim*dim
     * - else: just R
     */
    public double[] getRotations() { //TODO check caching
        if (valuesRinv == null) {
            return valuesR;
        }
        double[] out = new double[2 * dim * dim];
        System.arraycopy(valuesR, 0, out, 0, dim * dim);
        System.arraycopy(valuesRinv, 0, out, dim * dim, dim * dim);
        return out;
    }

    public double[] getValuesD() {
        return valuesD;
    }

    public BlockStructure getBlockStructure() {
        return blockStructure;
    }

    public boolean hasBlockStructure() {
        return blockStructure != null;
    }
}
