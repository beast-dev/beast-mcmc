package test.dr.math;

import dr.math.MathUtils;
import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

import static dr.math.ModifiedBesselFirstKind.*;

public class ModifiedBesselFirstKindTest extends TestCase {

    private static final double TOL = 1e-2;


    public void testModifiedBesselReal() {

        double[] xs = new double[]{1e-2, 1e-1, 1, 1e1, 1e2};
        int[] orders = new int[]{1, 2, 5, 10};

        double[] ys = new double[xs.length];
        for (int i = 0; i < xs.length; i++) {
            ys[i] = 1 / xs[i];
        }

        for (int i = 0; i < xs.length; i++) {
            for (int j = 0; j < orders.length; j++) {

                double bInt = bessi(xs[i], orders[j]);
                double bReal = bessi(xs[i], (double) orders[j]);

                if (!MathUtils.isRelativelyClose(bInt, bReal, TOL)) {
                    throw new AssertionFailedError("failed");
                }

                double byReal = bessi(ys[i], (double) orders[j]);
                double ratio = bessIRatio(xs[i], ys[i], orders[j]);

                if (!MathUtils.isRelativelyClose(bReal / byReal, ratio, TOL)) {
                    throw new AssertionFailedError("failed");
                }
            }
        }


    }


}
