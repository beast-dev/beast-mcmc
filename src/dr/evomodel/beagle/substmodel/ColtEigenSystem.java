package dr.evomodel.beagle.substmodel;

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

//        System.err.println("length = "+matrix.length + ","+matrix[0].length);

        EigenvalueDecomposition eigenDecomp = new EigenvalueDecomposition(new DenseDoubleMatrix2D(matrix));

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
//        DoubleMatrix1D eigenVImag = eigenDecomp.getImagEigenvalues();
        DoubleMatrix2D eigenVInv;

//        if (alegbra.cond(eigenV) > maxConditionNumber) {
//            return null;
//        }

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            return null;
        }

        double[][] Ievc = eigenVInv.toArray();
        double[][] Evec = eigenV.toArray();
        double[]   Eval = eigenVReal.toArray();
//        double[]   EvalImag = eigenVImag.toArray();

        return new EigenDecomposition(Evec,Ievc,Eval);
    }

    private static final double minProb = Property.DEFAULT.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);

}
