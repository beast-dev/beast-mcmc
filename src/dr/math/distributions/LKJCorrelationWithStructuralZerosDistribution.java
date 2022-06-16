package dr.math.distributions;

import dr.math.MathUtils;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;


public class LKJCorrelationWithStructuralZerosDistribution extends LKJCorrelationDistribution implements RandomGenerator {

    private final int[] blockAssignments;

    public static final String LKJ_WITH_ZEROS = "LKJCorrelationWithZerosDistribution";

    public static boolean DEBUG = false;


    public LKJCorrelationWithStructuralZerosDistribution(int dim, double shape, ArrayList<int[]> zeroBlocks) {
        super(dim, shape);

        this.blockAssignments = new int[dim];
        for (int i = 0; i < zeroBlocks.size(); i++) {
            for (int j = 0; j < zeroBlocks.get(i).length; j++) {
                blockAssignments[zeroBlocks.get(i)[j]] = i + 1; //want to save 0 for indices that aren't in any block
            }
        }

        System.err.println("Warning: LKJCorrelationDistribution with structural zeros does not have proper normalization constant"); //TODO
    }

    @Override
    public double[] nextRandom() {

        DenseMatrix64F partial = new DenseMatrix64F(dim, dim);

        for (int row = 0; row < dim; row++) {
            for (int col = row + 1; col < dim; col++) {

                if (blockAssignments[row] == 0 || blockAssignments[row] != blockAssignments[col]) {
                    int diag = col - row;
                    double alpha = shape + 0.5 * (dim - 1 - diag);
                    double beta = MathUtils.nextBeta(alpha, alpha);
                    beta *= 2;
                    beta -= 1;
                    partial.set(row, col, beta);
                    partial.set(col, row, beta);
                }

            }
        }

        if (DEBUG) {
            System.out.println(partial);
        }


        for (int i = 0; i < dim; i++) {
            partial.set(i, i, 1);
        }

        //convert partials to correlation matrix
        for (int diag = 2; diag < dim; diag++) {
            int dimSub = diag - 1;
            DenseMatrix64F Rinv = new DenseMatrix64F(dimSub, dimSub);
            DenseMatrix64F R = new DenseMatrix64F(dimSub, dimSub);

            DenseMatrix64F r1 = new DenseMatrix64F(dimSub, 1);
            DenseMatrix64F r2 = new DenseMatrix64F(dimSub, 1);

            DenseMatrix64F RInvr1 = new DenseMatrix64F(dimSub, 1);
            DenseMatrix64F RInvr2 = new DenseMatrix64F(dimSub, 1);

            for (int row = 0; row < dim - diag; row++) {
                int col = row + diag;

                for (int i = 0; i < dimSub; i++) {
                    int rowi = i + row + 1;
                    r1.set(i, 0, partial.get(rowi, row));
                    r2.set(i, 0, partial.get(rowi, col));
                    for (int j = 0; j < dimSub; j++) {
                        R.set(i, j, partial.get(rowi, j + row + 1));
                    }
                }

                CommonOps.invert(R, Rinv);
                CommonOps.mult(Rinv, r1, RInvr1);
                CommonOps.mult(Rinv, r2, RInvr2);

                double r1tRinvr1 = 0;
                double r1tRinvr2 = 0;
                double r2tRinvr2 = 0;

                for (int i = 0; i < dimSub; i++) {
                    r1tRinvr1 += RInvr1.get(i, 0) * r1.get(i, 0);
                    r1tRinvr2 += RInvr2.get(i, 0) * r1.get(i, 0);
                    r2tRinvr2 += RInvr2.get(i, 0) * r2.get(i, 0);
                }

                double d = (1 - r1tRinvr1) * (1 - r2tRinvr2);
                double c = r1tRinvr2 + partial.get(row, col) * Math.sqrt(d);
                partial.set(row, col, c);
                partial.set(col, row, c);
            }
        }

        double[] correlation = new double[upperTriangularSize(dim)];
        int ind = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = (i + 1); j < dim; j++) {
                correlation[ind] = partial.get(i, j);
                ind++;
            }
        }

        if (DEBUG) {
            System.out.println(partial);
        }

        return correlation;
    }

    @Override
    public double logPdf(Object x) {
        return logPdf((double[]) x);
    }

}
