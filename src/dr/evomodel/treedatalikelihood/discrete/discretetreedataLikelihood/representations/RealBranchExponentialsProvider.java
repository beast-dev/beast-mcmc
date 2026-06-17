package dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations;

/**
 * Optional zero-copy access to cached exp(t * lambda_i) values for
 * all-real spectral eigensystems.
 *
 * The borrowed slice is only valid until the next call that may mutate the
 * provider's branch-coefficient cache. Callers must consume it immediately.
 */
public interface RealBranchExponentialsProvider {

    boolean borrowRealBranchExponentials(int nodeNumber,
                                         double effectiveBranchLength,
                                         BorrowedSlice out);

    final class BorrowedSlice {
        private double[] values;
        private int offset;
        private int length;

        public double[] values() {
            return values;
        }

        public int offset() {
            return offset;
        }

        public int length() {
            return length;
        }

        void set(double[] values, int offset, int length) {
            this.values = values;
            this.offset = offset;
            this.length = length;
        }

        void clear() {
            this.values = null;
            this.offset = 0;
            this.length = 0;
        }
    }
}
