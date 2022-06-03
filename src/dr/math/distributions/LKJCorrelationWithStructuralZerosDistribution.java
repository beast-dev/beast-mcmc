package dr.math.distributions;

import dr.math.MathUtils;

import java.util.ArrayList;

public class LKJCorrelationWithStructuralZerosDistribution extends LKJCorrelationDistribution implements RandomGenerator {

    private final int[] blockAssignments;

    public static final String LKJ_WITH_ZEROS = "LKJCorrelationWithZerosDistribution";


    public LKJCorrelationWithStructuralZerosDistribution(int dim, double shape, ArrayList<int[]> zeroBlocks) {
        super(dim, shape);

        this.blockAssignments = new int[upperTriangularSize(dim)];
        for (int i = 0; i < zeroBlocks.size(); i++) {
            for (int j = 0; j < zeroBlocks.get(i).length; j++) {
                blockAssignments[zeroBlocks.get(i)[j]] = i + 1; //want to save 0 for indices that aren't in any block
            }
        }

        System.err.println("Warning: LKJCorrelationDistribution with structural zeros does not have proper normalization constant"); //TODO
    }

    @Override
    public double[] nextRandom() {
        int ind = 0;
        int n = upperTriangularSize(dim);
        double[] partialDraw = new double[n];

        for (int row = 0; row < dim; row++) {
            for (int col = row + 1; col < dim; col++) {

                if (blockAssignments[row] == 0 || blockAssignments[row] != blockAssignments[col]) {
                    int diag = row - col;
                    double alpha = shape + 0.5 * (dim - 1 - diag);
                    double beta = MathUtils.nextBeta(alpha, alpha);
                    beta *= 2;
                    beta -= 1;

                    partialDraw[ind] = beta;
                }

                ind++;
            }
        }

        double[] correlation = new double[n];

        //convert partials to correlation matrix
        for (int diag = 1; diag < dim; diag++) {
            //TODO
        }

        return correlation;
    }

    @Override
    public double logPdf(Object x) {
        return logPdf((double[]) x);
    }
}
