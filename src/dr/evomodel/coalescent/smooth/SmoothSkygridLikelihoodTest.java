package dr.evomodel.coalescent.smooth;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.TestCase;
import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

import java.util.ArrayList;
import java.util.List;

public class SmoothSkygridLikelihoodTest extends TestCase {

    private List<Tree> trees;
    private Parameter logPopSizeParameter;
    private Parameter gridPointParameter;
    private Parameter smoothRate;
    private GlobalSigmoidSmoothFunction smoothFunction = new GlobalSigmoidSmoothFunction();

    private final static UnivariateRealIntegrator integrator = new RombergIntegrator();

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        NewickImporter importer = new NewickImporter("(0:2.0,(1:1.0,2:1.0):1.0);");
        this.trees = new ArrayList<>();
        this.trees.add(importer.importTree(null));

        this.logPopSizeParameter = new Parameter.Default("logPopSize", new double[]{2.0, 1.5, 2.2});
        this.gridPointParameter = new Parameter.Default("logPopSize", new double[]{0.8, 1.6});

        this.smoothRate = new Parameter.Default(20.0);
    }

    public void testNewFormulation() throws Exception {
        SmoothSkygridLikelihood likelihood = new SmoothSkygridLikelihood("SmoothSkygridLikelihoodTest",
                trees, logPopSizeParameter, gridPointParameter, smoothRate);

        Tree tree = trees.get(0);
        final double startTime = 0.0;
        final double endTime = tree.getNodeHeight(tree.getRoot());

        final double stepLocation1 = tree.getNodeHeight(tree.getNode(1));
        final double stepLocation2 = tree.getNodeHeight(tree.getNode(3));

        GlobalSigmoidSmoothFunction sigmoidSmoothFunction = new GlobalSigmoidSmoothFunction();
        final double analytic = sigmoidSmoothFunction.getDoubleProductIntegration(startTime, endTime, stepLocation1,
                stepLocation2, smoothRate.getParameterValue(0));

        UnivariateRealFunction f = v -> getPairSigmoidProduct(v, stepLocation1, stepLocation2, smoothRate.getParameterValue(0));
        final double numeric = integrator.integrate(f, startTime, endTime);

        double diff = analytic - numeric;
    }

    public static double getPairSigmoidProduct(double time,
                                                 double stepLocation1, double stepLocation2,
                                                 double smoothRate) {
        GlobalSigmoidSmoothFunction sigmoidSmoothFunction = new GlobalSigmoidSmoothFunction();
        return sigmoidSmoothFunction.getSmoothValue(time, stepLocation1, 0, 1, smoothRate) *
                sigmoidSmoothFunction.getSmoothValue(time, stepLocation2, 0, 1, smoothRate);
    }

    public void testSmoothSkygridLikelihood () throws ConvergenceException, FunctionEvaluationException {
        SmoothSkygridLikelihood likelihood = new SmoothSkygridLikelihood("SmoothSkygridLikelihoodTest",
                trees, logPopSizeParameter, gridPointParameter, smoothRate);

        Tree tree = trees.get(0);
        double startTime = 0.0;
        double endTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double stepLocation1 = tree.getNodeHeight(tree.getNode(i));
            final double preStepValue1 = 0;
            final double postStepValue1 = i < tree.getExternalNodeCount() ? 1 : -1;
            for (int j = 0; j < tree.getNodeCount(); j++) {
                final double stepLocation2 = tree.getNodeHeight(tree.getNode(j));
                final double preStepValue2 = j == 0 ? -1 : 0;
                final double postStepValue2 = (j < tree.getExternalNodeCount() ? 1 : -1) + preStepValue2;
                for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                    final double stepLocation3 = gridPointParameter.getParameterValue(k);
                    final double preStepValue3 = k == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(0)) : 0;
                    final double postStepValue3 = k == 0? Math.exp(-logPopSizeParameter.getParameterValue(1)) :
                            Math.exp(-logPopSizeParameter.getParameterValue(k + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(k));
                    final double analytic = 0.5 * smoothFunction.getTripleProductIntegration(startTime, endTime,
                            stepLocation1, preStepValue1, postStepValue1,
                            stepLocation2, preStepValue2, postStepValue2,
                            stepLocation3, preStepValue3, postStepValue3,
                            smoothRate.getParameterValue(0));
                    sum += analytic;
                    UnivariateRealFunction f = v -> getTripleSigmoidProduct(v, stepLocation1, preStepValue1, postStepValue1, stepLocation2, preStepValue2, postStepValue2,stepLocation3, preStepValue3, postStepValue3, smoothRate.getParameterValue(0));
                    double numeric = 0.5 * integrator.integrate(f, 0.0, endTime);

                    assertEquals(analytic, numeric, 1e-6);
                }
            }
        }

        final double time = 0.6;
        double testValues[] = getSigmoidProductSum(time, tree);
        final double numericSum = likelihood.getReciprocalPopSizeInInterval(time, likelihood.getLineageCount(), likelihood.getPopulationSizeInverse());
        int nodeIndex = 4;
        final double smoothK = smoothFunction.getSmoothValue(time, tree.getNodeHeight(tree.getNode(nodeIndex)), 0, 1, smoothRate.getParameterValue(0));
        final double smoothKMinusOne = smoothFunction.getSmoothValue(time, tree.getNodeHeight(tree.getNode(nodeIndex)), -1, 0, smoothRate.getParameterValue(0));
        assertEquals(smoothK - 1, smoothKMinusOne, 1E-5);
        assertEquals(testValues[1], likelihood.getLineageCount().getLineageCount(time), 1e-5);
        assertEquals(testValues[2], likelihood.getLineageCount().getLineageCount(time) - 1, 1e-5);
        assertEquals(testValues[0], numericSum, 1e-5);


        final double lnL = likelihood.calculateLogLikelihood();
        UnivariateRealFunction f = v -> likelihood.getReciprocalPopSizeInInterval(v, likelihood.getLineageCount(), likelihood.getPopulationSizeInverse());
        final double integralBit = -integrator.integrate(f, 0.0, endTime);
        double multiplierBit = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(i + tree.getExternalNodeCount());
            multiplierBit += Math.log(likelihood.getPopulationSizeInverse().getPopulationSizeInverse(tree.getNodeHeight(node)));
        }
        assertEquals(lnL, integralBit + multiplierBit, 1e-6);

    }

    private double[] getSigmoidProductSum(double time, Tree tree) {
        double testSum = 0;
        double testK = 0;
        double kMinusOne = 0;

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double stepLocation1 = tree.getNodeHeight(tree.getNode(i));
            final double preStepValue1 = 0;
            final double postStepValue1 = i < tree.getExternalNodeCount() ? 1 : -1;
            testK += smoothFunction.getSmoothValue(time, stepLocation1, preStepValue1, postStepValue1,smoothRate.getParameterValue(0));
            kMinusOne = 0;
            for (int j = 0; j < tree.getNodeCount(); j++) {
                final double stepLocation2 = tree.getNodeHeight(tree.getNode(j));
                final double preStepValue2 = j == 0 ? -1 : 0;
                final double postStepValue2 = (j < tree.getExternalNodeCount() ? 1 : -1) + preStepValue2;
                kMinusOne += smoothFunction.getSmoothValue(time, stepLocation2, preStepValue2, postStepValue2, smoothRate.getParameterValue(0));
                for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                    final double stepLocation3 = gridPointParameter.getParameterValue(k);
                    final double preStepValue3 = k == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(0)) : 0;
                    final double postStepValue3 = k == 0? Math.exp(-logPopSizeParameter.getParameterValue(1)) :
                            Math.exp(-logPopSizeParameter.getParameterValue(k + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(k));
                    testSum += 0.5 * getTripleSigmoidProduct(time, stepLocation1, preStepValue1, postStepValue1, stepLocation2, preStepValue2, postStepValue2,
                            stepLocation3, preStepValue3, postStepValue3, smoothRate.getParameterValue(0));
                }
            }
        }
        return new double[]{testSum, testK, kMinusOne};
    }

    public static double getTripleSigmoidProduct(double time,
                                                 double stepLocation1, double preStepValue1, double postStepValue1,
                                                 double stepLocation2, double preStepValue2, double postStepValue2,
                                                 double stepLocation3, double preStepValue3, double postStepValue3,
                                                 double smoothRate) {
        GlobalSigmoidSmoothFunction sigmoidSmoothFunction = new GlobalSigmoidSmoothFunction();
        return sigmoidSmoothFunction.getSmoothValue(time, stepLocation1, preStepValue1, postStepValue1, smoothRate)
                * sigmoidSmoothFunction.getSmoothValue(time, stepLocation2, preStepValue2, postStepValue2, smoothRate)
                * sigmoidSmoothFunction.getSmoothValue(time, stepLocation3, preStepValue3, postStepValue3, smoothRate);
    }
}
