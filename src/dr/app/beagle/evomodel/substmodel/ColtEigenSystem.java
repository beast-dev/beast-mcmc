package dr.app.beagle.evomodel.substmodel;

import cern.colt.matrix.linalg.Property;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix2D;
import dr.math.matrixAlgebra.RobustEigenDecomposition;

/**
 * @author Marc Suchard
 */
public class ColtEigenSystem implements EigenSystem {

    public EigenDecomposition decomposeMatrix(double[][] matrix) {

        RobustEigenDecomposition eigenDecomp = new RobustEigenDecomposition(new DenseDoubleMatrix2D(matrix));

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix2D eigenVInv;

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            return null;
        }

        double[][] Evec = eigenV.toArray();
        double[][] Ievc = eigenVInv.toArray();
        double[] Eval = getAllEigenValues(eigenDecomp);

        double[] flatEvec = new double[Evec.length * Evec.length];
        double[] flatIevc = new double[Ievc.length * Ievc.length];

        for (int i = 0; i < Evec.length; i++) {
            System.arraycopy(Evec[i], 0, flatEvec, i * Evec.length, Evec.length);
            System.arraycopy(Ievc[i], 0, flatIevc, i * Ievc.length, Ievc.length);
        }

        return new EigenDecomposition(flatEvec, flatIevc, Eval);
    }

    protected double[] getAllEigenValues(RobustEigenDecomposition decomposition) {
        return decomposition.getRealEigenvalues().toArray();
    }

    private static final double minProb = Property.DEFAULT.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);

}
