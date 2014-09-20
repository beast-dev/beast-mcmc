package test.dr.evomodel.coalescent;

import junit.framework.TestCase;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evolution.util.Units;

/**
 * @author Joseph Heled
 *         Date: 23/06/2009
 */
public class VDdemographicFunctionTest extends TestCase {
    public void testExp() {
        // test that numerical and exact integration match (up to a point, numerical is not that good for those 
        // exponential gone to constant transitions.
        {
            double[] times = {1, 3};
            double[] logPops = {0, 2, 3};
            VariableDemographicModel.Type type = VariableDemographicModel.Type.EXPONENTIAL;
            Units.Type units = Units.Type.SUBSTITUTIONS;
            final VDdemographicFunction f = new VDdemographicFunction(times, logPops, units, type);

            double[][] vals = {{0, 2, 1e-9}, {1, 2, 1e-9}, {0, 5, 1e-6}, {2, 6, 1e-6}};
            for(double[] c : vals) {
                double v1 = f.getIntegral(c[0], c[1]);
                double v2 = f.getNumericalIntegral(c[0], c[1]);
                assertEquals(Math.abs(1 - v1 / v2), 0, c[2]);
            }
        }

        {
            double[] times = {1, 3};
            // try a const interval
            double[] logPops = {0, 0, 3};
            VariableDemographicModel.Type type = VariableDemographicModel.Type.EXPONENTIAL;
            Units.Type units = Units.Type.SUBSTITUTIONS;
            final VDdemographicFunction f = new VDdemographicFunction(times, logPops, units, type);

            double[][] vals = {{0, .7, 1e-9}, {1, 2, 1e-9}, {0, 5, 1e-6}, {2, 6, 1e-6}};
            for(double[] c : vals) {
                double v1 = f.getIntegral(c[0], c[1]);
                double v2 = f.getNumericalIntegral(c[0], c[1]);
                assertEquals(Math.abs(1 - v1 / v2), 0, c[2]);
            }
        }
    }
}
