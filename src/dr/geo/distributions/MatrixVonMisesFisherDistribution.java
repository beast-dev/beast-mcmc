package dr.geo.distributions;

import dr.math.MathUtils;
import dr.math.distributions.RandomGenerator;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;

import static dr.math.MathUtils.nextDouble;

public class MatrixVonMisesFisherDistribution implements RandomGenerator {

    private final DenseMatrix64F C;
    private final int nRows;
    private final int nColumns;
    private final DenseMatrix64F mkBuffer1;
    private final DenseMatrix64F mkBuffer2;
    private final DenseMatrix64F kkBuffer1;
    private final DenseMatrix64F kkBuffer2;
    private final DenseMatrix64F kkBuffer3;
    private static final int MAX_REJECTS = 100;

    public MatrixVonMisesFisherDistribution(int nRows, int nColumns) {
        this.C = new DenseMatrix64F(nRows, nColumns);
        this.nRows = nRows;
        this.nColumns = nColumns;
        this.mkBuffer1 = new DenseMatrix64F(nRows, nColumns);
        this.mkBuffer2 = new DenseMatrix64F(nRows, nColumns);
        this.kkBuffer1 = new DenseMatrix64F(nColumns, nColumns);
        this.kkBuffer2 = new DenseMatrix64F(nColumns, nColumns);
        this.kkBuffer3 = new DenseMatrix64F(nColumns, nColumns);
    }


    @Override
    public double[] nextRandom() {

        int rejects = 0;

        while (rejects < MAX_REJECTS) {
            DenseMatrix64F uniformDraw = nextUniform();

            SingularValueDecomposition svd = new DecompositionFactory().svd(C.numRows, C.numCols, false, false, true);
            svd.decompose(C);
            double[] singularValues = svd.getSingularValues();
            CommonOps.multTransA(C, uniformDraw, kkBuffer1); //TODO: just need the trace, super inefficient
            double trace = 0;
            for (int i = 0; i < C.numCols; i++) {
                trace += kkBuffer1.get(i, i) - singularValues[i];
            }

            if (nextDouble() < trace) {
                return uniformDraw.getData();
            }
        }

        throw new RuntimeException("Rejection sampler failed.");//TODO: handle better

    }

    private DenseMatrix64F nextUniform() {
        double[] X = new double[C.getNumElements()];
        for (int i = 0; i < X.length; i++) {
            X[i] = MathUtils.nextGaussian();
        }
        mkBuffer1.setData(X);
        CommonOps.multTransB(mkBuffer1, mkBuffer1, kkBuffer1);
        EigenDecomposition<DenseMatrix64F> eig = DecompositionFactory.eig(nColumns, true, true);
        eig.decompose(kkBuffer1);
        for (int i = 0; i < nColumns; i++) {
            DenseMatrix64F vector = eig.getEigenVector(i);
            for (int j = 0; j < nColumns; j++) {
                kkBuffer2.set(j, i, vector.get(j, 0));
            }
            double value = eig.getEigenvalue(i).getReal(); //TODO: should be real
            kkBuffer3.set(i, i, 1.0 / Math.sqrt(value));
        }

        CommonOps.mult(kkBuffer2, kkBuffer3, kkBuffer1); // kkBuffer1 = VD
        CommonOps.multTransB(kkBuffer1, kkBuffer2, kkBuffer3); // kkBuffer3 = VDVt
        CommonOps.mult(mkBuffer1, kkBuffer3, mkBuffer2);
        return mkBuffer2;
    }


    @Override
    public double logPdf(Object x) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void setC(double[] values) {
        C.setData(values);
    }
}
