package test.dr.inference.model;

import dr.inference.model.Parameter;
import dr.inference.model.VectorSliceParameter;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */
public class VectorSliceParameterTest extends MathTestCase {
    public void testVectorSlice() {

        Parameter pA = new Parameter.Default(new double[] { 0.0, 1.0, 2.0 });
        Parameter pB = new Parameter.Default(new double[] { 3.0, 4.0, 5.0 });

        VectorSliceParameter slice = new VectorSliceParameter("slice",1);
        slice.addParameter(pA);
        slice.addParameter(pB);

        assertEquals(2, slice.getDimension());
        assertEquals(1.0, slice.getParameterValue(0), tolerance);
        assertEquals(4.0, slice.getParameterValue(1), tolerance);

        pA.setParameterValue(1, 41.0);
        assertEquals(41.0, slice.getParameterValue(0), tolerance);

        slice.setParameterValue(1, 44.0);
        assertEquals(44.0, pB.getParameterValue(1), tolerance);

    }
    private static final double tolerance = 1E-6;
}
