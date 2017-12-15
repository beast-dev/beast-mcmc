package dr.evomodel.treedatalikelihood.continuous.cdi;

import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @author Paul Bastide
 */

public class SafeMultivariateActualizedWithDriftIntegratorTest {
//    @Test
//    public void getBranchMatrices() throws Exception {
//    }

    interface Instance {
        int getDimTrait();

        double[] getPrecision();

        double[] getSelectionStrength();

        double[] getStationaryVariance();
    }

    /**
     *     Results obtained using the following R code:
     *
     *     variance <- matrix(data = c(1, 0, 0, 1), 2, 2)
     *     alpha <- matrix(data = c(1, 0, 0, 1), 2, 2)
     *     gamma <- PylogeneticEM::compute_stationary_variance(variance, alpha)
     *
     */

    Instance test0 = new Instance() {
        public int getDimTrait(){ return 2;}

        public double[] getPrecision() {
            return new double[]{
                    1.0, 0.0,
                    0.0, 1.0
            };
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    1.0, 0.0,
                    0.0, 1.0
            };
        }

        public double[] getStationaryVariance() {
            return new double[]{
                    0.5, 0.0,
                    0.0, 0.5
            };
        }
    };

    Instance test1 = new Instance() {
        public int getDimTrait(){ return 2;}

        public double[] getPrecision() {
            return new double[]{
                    1.4, 0.4,
                    0.4, 1.0
            };
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    2.3, 0.2,
                    0.2, 1.2
            };
        }

        public double[] getStationaryVariance() {
            return new double[]{
                    0.1867037, -0.1309637,
                    -0.1309637,  0.4922574
            };
        }
    };

    Instance test2 = new Instance() {
        public int getDimTrait(){ return 5;}

        public double[] getPrecision() {
            return new double[]{
                    1.4, 0.4, 0.0, 0.1, 0.9,
                    0.4, 1.0, 0.3, 0.0, 0.0,
                    0.0, 0.3, 2.1, 0.2, 0.2,
                    0.1, 0.0, 0.2, 1.2, 0.5,
                    0.9, 0.0, 0.2, 0.5, 2.0
            };
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    1.4, 0.0, 0.0, 0.0, 1.0,
                    0.0, 2.0, 0.5, 0.0, 0.0,
                    0.0, 0.5, 3.1, 0.0, 0.0,
                    0.0, 0.0, 0.0, 2.2, 0.4,
                    1.0, 0.0, 0.0, 0.4, 1.0
            };
        }

        public double[] getStationaryVariance() {
            return new double[]{
                     2.55009871, -0.21879549,  0.064601442,  0.578712147, -2.94641261,
                    -0.21879549,  0.33965614, -0.082974248, -0.024285677,  0.17602648,
                     0.06460144, -0.08297425,  0.097219008, -0.004240138, -0.05931872,
                     0.57871215, -0.02428568, -0.004240138,  0.357880727, -0.77676809,
                    -2.94641261,  0.17602648, -0.059318715, -0.776768090,  3.68424701
            };
        }
    };

    Instance[] all = {test0, test1, test2};

    @Test
    public void setDiffusionStationaryVariance() throws Exception {
        for (Instance test : all) {
            PrecisionType precisionType = PrecisionType.FULL;
            int numTraits = 1;
            int dimTrait = test.getDimTrait();
            int partialBufferCount = 1;
            int matrixBufferCount = 1;

            SafeMultivariateActualizedWithDriftIntegrator cdi = new SafeMultivariateActualizedWithDriftIntegrator(
                    precisionType,
                    numTraits,
                    dimTrait,
                    partialBufferCount,
                    matrixBufferCount
            );

            double precision[] = test.getPrecision();
            cdi.setDiffusionPrecision(0, precision, 2.0);

            double alpha[] = test.getSelectionStrength();

            cdi.setDiffusionStationaryVariance(0, alpha);

            double[] stationaryVariance = cdi.getStationaryVariance(0);

            double expectedStationaryVariance[] = test.getStationaryVariance();

            assertArrayEquals(stationaryVariance, expectedStationaryVariance, 1e-6);
        }
    }


//    @Test
//    public void updateOrnsteinUhlenbeckDiffusionMatrices() throws Exception {
//    }
//
//    @Test
//    public void updatePreOrderPartial() throws Exception {
//    }
//
//    @Test
//    public void updatePartial() throws Exception {
//    }

}