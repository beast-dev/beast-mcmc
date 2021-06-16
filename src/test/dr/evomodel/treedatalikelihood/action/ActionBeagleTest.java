package test.dr.evomodel.treedatalikelihood.action;

import dr.evomodel.treedatalikelihood.action.ActionBeagleDelegate;
import org.newejml.data.DMatrixSparseCSC;
import org.newejml.data.DMatrixSparseTriplet;
import org.newejml.ops.DConvertMatrixStruct;
import test.dr.math.MathTestCase;

public class ActionBeagleTest extends MathTestCase {
    private ActionBeagleDelegate beagle;
    private double[] stationaryFrequency;
    private DMatrixSparseTriplet Q;

    private int patternCount = 2;
    private int tipCount = 3;
    private int stateCount = 4;
    private int categoryCount = 1;
    private int matrixBufferCount = 4;
    private int partialsSize;

    public void setUp() throws Exception {

        this.partialsSize = patternCount * stateCount * categoryCount;

        this.stationaryFrequency = new double[]{0.1, 0.3, 0.2, 0.4};

        this.Q = new DMatrixSparseTriplet(4, 4, 16);

        for (int i = 0; i < 4; i++) {
            double sum = 0.0;
            for (int j = 0; j < 4; j++) {
                if (i == j) continue;
                Q.addItem(i, j, stationaryFrequency[j]);
                sum += stationaryFrequency[j];
            }
            Q.addItem(i, i, -sum);
        }
    }

    public void testQsetup() {
        DMatrixSparseCSC Qc = DConvertMatrixStruct.convert(Q, (DMatrixSparseCSC) null);
        this.beagle = new ActionBeagleDelegate(tipCount, (2 * tipCount - 2) * 2, patternCount,
                stateCount, categoryCount, matrixBufferCount, partialsSize, new DMatrixSparseCSC[]{Qc.copy(), Qc.copy(), Qc.copy(), Qc.copy()});

        beagle.setPartials(0, new double[]{
                1., 0., 0., 0.,
                1., 0., 0., 0.,
                1., 0., 0., 0.,
                0., 0., 0., 1.,
                0., 1., 0., 0.});
        beagle.setPartials(1, new double[]{
                0., 0., 1., 0.,
                1., 0., 0., 0.,
                0., 0., 1., 0.,
                0., 0., 0., 1.,
                1., 0., 0., 0.
        });
        beagle.setPartials(2, new double[]{
                0., 0., 1., 0.,
                1., 0., 0., 0.,
                0., 0., 1., 0.,
                0., 0., 1., 0.,
                1., 0., 0., 0.
        });
        int[] operations = new int[]{
                3, 0, 0, 0, 0, 1, 1
        };
        beagle.updatePartials(operations, 1, 0);



    }


}
