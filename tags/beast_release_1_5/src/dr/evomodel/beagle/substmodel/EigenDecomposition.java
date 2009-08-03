package dr.evomodel.beagle.substmodel;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @Author Marc A. Suchard
 * @version $Id$
 */
public class EigenDecomposition {

    public EigenDecomposition(double[] evec, double[] ievc, double[] eval) {
        Evec = evec;
        Ievc = ievc;
        Eval = eval;
    }

    public EigenDecomposition copy() {
        double[] evec = Evec.clone();
        double[] ievc = Ievc.clone();
        double[] eval = Eval.clone();

        return new EigenDecomposition(evec, ievc, eval);
    }

    /**
     * This function returns the Eigen vectors.
     * @return the array
     */
    public final double[] getEigenVectors() {
        return Evec;
    }

    /**
     * This function returns the inverse Eigen vectors.
     * @return the array
     */
    public final double[] getInverseEigenVectors() {
        return Ievc;
    }

    /**
     * This function returns the Eigen values.
     * @return the Eigen values
     */
    public final double[] getEigenValues() {
        return Eval;
    }

    // Eigenvalues, eigenvectors, and inverse eigenvectors
    private final double[] Evec;
    private final double[] Ievc;
    private final double[] Eval;

}
