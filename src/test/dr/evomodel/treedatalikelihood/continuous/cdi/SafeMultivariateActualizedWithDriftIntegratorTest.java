package test.dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.EigenOps;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

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

        double[] getEigenValuesStrengthOfSelection();

        double[] getEigenVectorsStrengthOfSelection();

        boolean isASymmetric();

        abstract class Basic implements Instance {

            private EigenDecomposition decomposeStrenghtOfSelection(double[] Aparam) {
                int n = getDimTrait();
                DenseMatrix64F A = MissingOps.wrap(Aparam, 0, n, n);
                // Decomposition
                EigenDecomposition eigA = DecompositionFactory.eig(n, true, false);
                if (!eigA.decompose(A)) throw new RuntimeException("Eigen decomposition failed.");
                return eigA;
            }

            final EigenDecomposition eigenDecompositionStrengthOfSelection = decomposeStrenghtOfSelection(getSelectionStrength());

            public double[] getEigenValuesStrengthOfSelection() {
                int dim = getDimTrait();
                double[] eigA = new double[dim];
                for (int p = 0; p < dim; ++p) {
                    Complex64F ev = eigenDecompositionStrengthOfSelection.getEigenvalue(p);
                    if (!ev.isReal())
                        throw new RuntimeException("Selection strength A should only have real eigenvalues.");
                    eigA[p] = ev.real;
                }
                return eigA;
            }

            public double[] getEigenVectorsStrengthOfSelection() {
                return EigenOps.createMatrixV(eigenDecompositionStrengthOfSelection).data;
            }

        }
    }

    /**
     * Results obtained using the following R code:
     * <p>
     * variance <- solve(matrix(data = c(1, 0, 0, 1), 2, 2))
     * alpha <- matrix(data = c(1, 0, 0, 1), 2, 2)
     * gamma <- PhylogeneticEM::compute_stationary_variance(variance, alpha)
     */

    Instance test0 = new Instance.Basic() {
        public int getDimTrait() {
            return 2;
        }

        public boolean isASymmetric() {
            return true;
        }

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

    Instance test1 = new Instance.Basic() {
        public int getDimTrait() {
            return 2;
        }

        public boolean isASymmetric() {
            return true;
        }

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
                    -0.1309637, 0.4922574
            };
        }
    };

    Instance test2 = new Instance.Basic() {
        public int getDimTrait() {
            return 5;
        }

        public boolean isASymmetric() {
            return true;
        }

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
                    2.55009871, -0.21879549, 0.064601442, 0.578712147, -2.94641261,
                    -0.21879549, 0.33965614, -0.082974248, -0.024285677, 0.17602648,
                    0.06460144, -0.08297425, 0.097219008, -0.004240138, -0.05931872,
                    0.57871215, -0.02428568, -0.004240138, 0.357880727, -0.77676809,
                    -2.94641261, 0.17602648, -0.059318715, -0.776768090, 3.68424701
            };
        }
    };

    Instance test3 = new Instance.Basic() {
        public int getDimTrait() {
            return 2;
        }

        public boolean isASymmetric() {
            return false;
        }

        public double[] getPrecision() {
            return new double[]{
                    1.4, 0.4,
                    0.4, 1.0
            };
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    2.3, 1.0,
                    0.2, 1.2
            };
        }

        public double[] getStationaryVariance() {
            return new double[]{
                    0.2862183, -0.2550763,
                    -0.2550763,  0.5129428
            };
        }
    };

    Instance[] all = {test0, test1, test2, test3};

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
                    dimTrait,
                    partialBufferCount,
                    matrixBufferCount,
                    test.isASymmetric()
            );

            double precision[] = test.getPrecision();
            cdi.setDiffusionPrecision(0, precision, 2.0);

            double[] alphaEig = test.getEigenValuesStrengthOfSelection();
            double[] alphaRot = test.getEigenVectorsStrengthOfSelection();

            cdi.setDiffusionStationaryVariance(0, alphaEig, alphaRot);

            double[] stationaryVariance = cdi.getStationaryVariance(0);

            double expectedStationaryVariance[] = test.getStationaryVariance();

            assertEquals(stationaryVariance.length, expectedStationaryVariance.length, 1e-6);
            for (int i = 0; i < expectedStationaryVariance.length; ++i) {
                assertEquals(stationaryVariance[i], expectedStationaryVariance[i], 1e-6);
            }
        }
    }


//    @Test
//    public void updateOrnsteinUhlenbeckDiffusionMatrices() throws Exception {
//        for (Instance test : all) {
//            PrecisionType precisionType = PrecisionType.FULL;
//            int numTraits = 1;
//            int dimTrait = test.getDimTrait();
//            int partialBufferCount = 1;
//            int matrixBufferCount = 1;
//
//            SafeMultivariateActualizedWithDriftIntegrator cdi = new SafeMultivariateActualizedWithDriftIntegrator(
//                    precisionType,
//                    numTraits,
//                    dimTrait,
//                    partialBufferCount,
//                    matrixBufferCount
//            );
//
//            double precision[] = test.getPrecision();
//            cdi.setDiffusionPrecision(0, precision, 2.0);
//
//            double alpha[] = test.getSelectionStrength();
//
//            cdi.setDiffusionStationaryVariance(0, alpha);
//
//            final double[] edgeLengths = test.getEdgeLengths();
//            final double[] optimalRates = test.getOptimalRates();
//
//            cdi.updateOrnsteinUhlenbeckDiffusionMatrices(0, 0,
//                    edgeLengths,
//                    optimalRates,
//                    alpha,
//                    1)
//
//        }
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