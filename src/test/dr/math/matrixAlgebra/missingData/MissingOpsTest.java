package test.dr.math.matrixAlgebra.missingData;

import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Paul Bastide
 */

public class MissingOpsTest {

    interface Instance {

        int getRank();

        double getLogDet();

        DenseMatrix64F getMatrix();

        abstract class Basic implements Instance {

        }
    }

    Instance test0 = new Instance.Basic() {
        public int getRank() {
            return 2;
        }

        public double getLogDet() {
            return 42.160726233513714;
        }

        public DenseMatrix64F getMatrix() {
            return new DenseMatrix64F(2, 2, true,
                    0.004255873897787696, 0.010962329615067505,
                           0.010962329615067505, 0.02823689645782389);
        }
    };

    Instance test1 = new Instance.Basic() {
        public int getRank() {
            return 2;
        }

        public double getLogDet() {
            return 42.56618093096408;
        }

        public DenseMatrix64F getMatrix() {
            return new DenseMatrix64F(2, 2, true,
                    0.004255918204391931, 0.010962443740752292,
                           0.010962443740752292, 0.028237190424652364);
        }
    };


    Instance[] all = {test0, test1};

    @Test
    public void safeDeterminant() throws Exception {
        for (Instance test : all) {

            int rank_test = test.getRank();
            double logDet_test = test.getLogDet();

            DenseMatrix64F P = test.getMatrix();

            InversionResult c = MissingOps.safeDeterminant(P, false);

            assertEquals(c.getLogDeterminant(), logDet_test, 1e-6);
            assertEquals(c.getEffectiveDimension(), rank_test, 1e-6);
        }
    }
}