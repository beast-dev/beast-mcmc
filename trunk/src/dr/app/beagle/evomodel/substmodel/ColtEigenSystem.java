package dr.app.beagle.evomodel.substmodel;

import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Property;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix1D;

/**
 * @author Marc Suchard
 */
public class ColtEigenSystem implements EigenSystem {

    public EigenDecomposition decomposeMatrix(double[][] matrix) {

        EigenvalueDecomposition eigenDecomp = new EigenvalueDecomposition(new DenseDoubleMatrix2D(matrix));

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
        DoubleMatrix2D eigenVInv;

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            return null;
        }

        double[][] Evec = eigenV.toArray();
        double[][] Ievc = eigenVInv.toArray();
        double[] Eval = eigenVReal.toArray();

        double[] flatEvec = new double[Evec.length * Evec.length];
        double[] flatIevc = new double[Ievc.length * Ievc.length];

        for (int i = 0; i < Evec.length; i++) {
            System.arraycopy(Evec[i], 0, flatEvec, i * Evec.length, Evec.length);
            System.arraycopy(Ievc[i], 0, flatIevc, i * Ievc.length, Ievc.length);
        }

        return new EigenDecomposition(flatEvec, flatIevc, Eval);
    }

    private static final double minProb = Property.DEFAULT.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);

}
