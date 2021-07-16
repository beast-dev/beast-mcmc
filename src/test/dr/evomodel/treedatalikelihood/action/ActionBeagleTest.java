package test.dr.evomodel.treedatalikelihood.action;

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.treedatalikelihood.action.ActionBeagleDelegate;
import dr.evomodel.treedatalikelihood.action.ActionEvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.action.HomogeneousActionSubstitutionModelDelegate;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import org.newejml.data.DMatrixSparseCSC;
import org.newejml.data.DMatrixSparseTriplet;
import org.newejml.ops.DConvertMatrixStruct;
import test.dr.math.MathTestCase;

public class ActionBeagleTest extends MathTestCase {
    private ActionBeagleDelegate beagle;
    private double[] stationaryFrequency;
    private DMatrixSparseTriplet Q;

    private int patternCount = 5;
    private int tipCount = 3;
    private int stateCount = 4;
    private int categoryCount = 2;
    private int matrixBufferCount = 4;
    private int partialsSize;
    private int partialsBufferCount;

    public void setUp() throws Exception {

        this.partialsSize = stateCount * patternCount * categoryCount;
        this.partialsBufferCount = (2 * tipCount - 2) * 2 ;


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
        Parameter kappa = new Parameter.Default(3.0);
        Parameter pi = new Parameter.Default(new double[]{0.1, 0.3, 0.2, 0.4});
        FrequencyModel frequencyModel = new FrequencyModel(Nucleotides.INSTANCE, pi);
        HKY hky = new HKY(kappa, frequencyModel);
        PartialsRescalingScheme rescalingScheme = PartialsRescalingScheme.AUTO;
        ActionEvolutionaryProcessDelegate evolutionaryProcessDelegate = new HomogeneousActionSubstitutionModelDelegate(hky, 5);
        this.beagle = new ActionBeagleDelegate(tipCount, partialsBufferCount, patternCount,
                stateCount, categoryCount, matrixBufferCount, partialsSize,
                rescalingScheme, evolutionaryProcessDelegate);
        evolutionaryProcessDelegate.updateSubstitutionModels(beagle, false);

        // nCategory = 1;
//        beagle.setPartials(0, new double[]{
//                1., 0., 0., 0.,
//                1., 0., 0., 0.,
//                1., 0., 0., 0.,
//                0., 0., 0., 1.,
//                0., 1., 0., 0.});
//        beagle.setPartials(1, new double[]{
//                0., 0., 1., 0.,
//                1., 0., 0., 0.,
//                0., 0., 1., 0.,
//                0., 0., 0., 1.,
//                1., 0., 0., 0.
//        });
//        beagle.setPartials(2, new double[]{
//                0., 0., 1., 0.,
//                1., 0., 0., 0.,
//                0., 0., 1., 0.,
//                0., 0., 1., 0.,
//                1., 0., 0., 0.
//        });
        // nCategory = 2;
        beagle.setPartials(0, new double[]{
                1., 0., 0., 0.,
                1., 0., 0., 0.,
                1., 0., 0., 0.,
                0., 0., 0., 1.,
                0., 1., 0., 0.,
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
                1., 0., 0., 0.,
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
                1., 0., 0., 0.,
                0., 0., 1., 0.,
                1., 0., 0., 0.,
                0., 0., 1., 0.,
                0., 0., 1., 0.,
                1., 0., 0., 0.
        });
        double[] patternWeights = new double[]{1.0, 1.0, 1.0, 1.0, 1.0};
        beagle.setPatternWeights(patternWeights);
        double[] categoryWeights = new double[]{0.2, 0.8};
        beagle.setCategoryWeights(0, categoryWeights);
        double[] categoryRates = new double[]{1.2, 0.5};
        beagle.setCategoryRates(categoryRates);
        evolutionaryProcessDelegate.updateTransitionMatrices(beagle,
                new int[]{0, 1, 2, 3},
                new double[]{0.1, 0.1, 0.2, 0.3},
                4,
                false);
        int[] operations = new int[]{
                3, 0, 0, 0, 0, 1, 1,
                4, 0, 0, 3, 3, 2, 2
        };
        beagle.updatePartials(operations, 2, 0);

        double[] seePartials = new double[partialsSize];
        beagle.getPartials(4, 0, seePartials);

        beagle.setStateFrequencies(0, stationaryFrequency);
        int rootIndex = 4;
        double[] sumLogLikelihoods = new double[1];
        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                new int[]{0}, 1, sumLogLikelihoods);

        double logL = sumLogLikelihoods[0];




    }


}
