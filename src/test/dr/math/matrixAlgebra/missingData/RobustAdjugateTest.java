package test.dr.math.matrixAlgebra.missingData;

import dr.math.matrixAlgebra.EJMLUtils;
import org.ejml.data.DenseMatrix64F;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RobustAdjugateTest {

    interface Instance {
        DenseMatrix64F getMatrix();

        double[] getTrueAdjugate();

        abstract class Basic implements Instance {

        }
    }

    Instance test0 = new Instance.Basic() {

        @Override
        public DenseMatrix64F getMatrix() {
            int dim = 4;
            return DenseMatrix64F.wrap(
                    dim, dim,
                    new double[]{0.9866964953013155, 0.15545690836486598, -1.1513262141217357, 0.4725542031241179,
                            -0.19421901601084784, 0.5978801523926415, -0.45827856969399533, -0.09859711476482771,
                            0.8510272497600995, 0.46551635730932206, -0.24277738841539084, 0.5870842918793533,
                            1.3516144337818743, 0.3483227354656835, 0.3519407851667098, -0.26191847550729874}
            );
        }

        @Override
        public double[] getTrueAdjugate() {
            return new double[]{-0.25959911412615355, 0.2806715010354877, 0.05816932938829323, -0.44364091890937585,
                    0.49990493839334016, -1.0147493713240305, -0.6250795751738031, -0.1171771471159761,
                    0.7038993830503742, 0.3213216171925368, -0.5630922852614357, -0.11313981508126116,
                    0.2710051185801451, 0.5306443430372886, -1.28773681446716, 0.689220622236884};
        }
    };

    Instance test1 = new Instance.Basic() {

        @Override
        public DenseMatrix64F getMatrix() {
            int dim = 4;
            return DenseMatrix64F.wrap(dim, dim,
                    new double[]{1, 0.0, -0.4564755436471619, 0.20065293735042314,
                            0.0, 1, -0.8800002808267111, 0.041910304168150225,
                            -0.4564755436471619, -0.8800002808267111, 1, 0.0,
                            0.20065293735042314, 0.041910304168150225, 0.0, 1.0});
        }

        @Override
        public double[] getTrueAdjugate() {
            return new double[]{0.22384303214944315, 0.4101080322366171, 0.46307405332955914, -0.06210251427904914,
                    0.4101080322366171, 0.7513684767846751, 0.8484087575441298, -0.11377946270368472,
                    0.46307405332955914, 0.8484087575441298, 0.9579819251371803, -0.12847423809893943,
                    -0.06210251427904914, -0.11377946270368472, -0.1284742380989394, 0.017229583796937953};
        }
    };

    private static final double TOL = 1e-12;


    @Test
    public void adjugate() throws Exception {

        Instance[] tests = new Instance[]{test0, test1};
        for (Instance test : tests) {
            DenseMatrix64F X = test.getMatrix();

            DenseMatrix64F A = EJMLUtils.computeRobustAdjugate(X);

            double[] a = test.getTrueAdjugate();

            for (int i = 0; i < a.length; i++) {
                assertEquals(A.get(i), a[i], TOL);
            }

        }
    }
}
