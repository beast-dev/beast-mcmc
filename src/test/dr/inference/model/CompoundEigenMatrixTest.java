package test.dr.inference.model;

import dr.inference.model.CompoundEigenMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
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

public class CompoundEigenMatrixTest {

    interface Instance {
        int getDimTrait();

        double[] getSelectionStrength();

        Parameter getEigenValuesStrengthOfSelection();

        MatrixParameter getEigenVectorsStrengthOfSelection();

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

            public Parameter getEigenValuesStrengthOfSelection() {
                int dim = getDimTrait();
                double[] eigA = new double[dim];
                for (int p = 0; p < dim; ++p) {
                    Complex64F ev = eigenDecompositionStrengthOfSelection.getEigenvalue(p);
                    if (!ev.isReal())
                        throw new RuntimeException("Selection strength A should only have real eigenvalues.");
                    eigA[p] = ev.real;
                }
                return new Parameter.Default(eigA);
            }

            public MatrixParameter getEigenVectorsStrengthOfSelection() {
                int dim = getDimTrait();
                DenseMatrix64F V = EigenOps.createMatrixV(eigenDecompositionStrengthOfSelection);

                Parameter[] eigenVectors = new Parameter[getDimTrait()];
                double[] column = new double[dim - 1];
                double norm;
                double sign;
                for (int i = 0; i < dim; i++) {
                    norm = 0.0;
                    sign = 1.0;
                    for (int j = 0; j < dim; j++) {
                        norm += V.get(j, i) * V.get(j, i);
                    }
                    norm = Math.sqrt(norm);
                    if (V.get(dim - 1, i) < 0) sign = -1.0;
                    for (int j = 0; j < dim - 1; j++) {
                        column[j] = V.get(j, i) / norm * sign;
                    }
                    eigenVectors[i] = new Parameter.Default(column);
                }
                return new MatrixParameter("attenuationMatrix", eigenVectors);
            }

        }
    }

    Instance test0 = new Instance.Basic() {
        public int getDimTrait() {
            return 2;
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    1.0, 0.0,
                    0.0, 1.0
            };
        }
    };

    Instance test1 = new Instance.Basic() {
        public int getDimTrait() {
            return 2;
        }


        public double[] getSelectionStrength() {
            return new double[]{
                    2.3, 0.2,
                    0.2, 1.2
            };
        }
    };

    Instance test2 = new Instance.Basic() {
        public int getDimTrait() {
            return 5;
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    1.4, 0.0, 0.0, 0.0, 1.0,
                    0.0, 2.0, 0.5, 0.0, 0.0,
                    0.3, 0.5, 3.1, 0.0, 0.0,
                    0.4, 0.0, 0.0, 2.2, 0.4,
                    1.0, 0.0, 0.0, 0.4, 1.0
            };
        }
    };

    Instance test3 = new Instance.Basic() {
        public int getDimTrait() {
            return 2;
        }

        public double[] getSelectionStrength() {
            return new double[]{
                    2.3, 1.0,
                    0.2, 1.2
            };
        }
    };

    Instance test4 = new Instance.Basic() {
        public int getDimTrait() {
            return 4;
        }

        public double[] getSelectionStrength() {
            return new double[]{
                      5.635055,  0.84322690,  2.200631,  1.6066434,
                     -2.123113, -0.07568064, -0.879367, -0.7606432,
                     -3.900102, -0.54334995, -1.397044, -1.1259163,
                    -11.481042, -1.89373481, -4.613851, -3.1623302
            };
        }
    };

    Instance[] all = {test0, test1, test2, test3, test4};

    @Test
    public void CompoundEigenMatrix() throws Exception {
        for (Instance test : all) {

            Parameter alphaEig = test.getEigenValuesStrengthOfSelection();
            MatrixParameter alphaRot = test.getEigenVectorsStrengthOfSelection();

            CompoundEigenMatrix alphaComp = new CompoundEigenMatrix(alphaEig, alphaRot);

            double[] alphaCompo = alphaComp.getParameterValues();

            double[] alphaCompoAttibute = alphaComp.getAttributeValue();

            double[] alphaExpected = test.getSelectionStrength();

            assertEquals(alphaCompo.length, alphaExpected.length, 1e-6);
            for (int i = 0; i < alphaCompo.length; ++i) {
                assertEquals("alpha " + i,
                        alphaCompo[i], alphaExpected[i], 1e-6);
            }

            for (int i = 0; i < alphaCompo.length; ++i) {
                assertEquals("alpha " + i,
                        alphaCompo[i], alphaCompoAttibute[i], 1e-6);
            }
        }
    }
}