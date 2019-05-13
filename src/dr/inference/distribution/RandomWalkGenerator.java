package dr.inference.distribution;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.demographicmodel.PiecewisePopulationModel;
import dr.evomodel.coalescent.TreeIntervals;
import dr.inference.model.*;
import dr.inferencexml.distribution.RandomWalkGeneratorParser;
import dr.math.MathUtils;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

/**
 * Created by mkarcher on 4/3/17.
 */
public class RandomWalkGenerator extends AbstractModelLikelihood implements GaussianProcessRandomGenerator {

    public RandomWalkGenerator(Parameter data, Parameter firstElementPrecision, Parameter precision) {
        super(RandomWalkGeneratorParser.RANDOM_WALK_GENERATOR);
        this.data = data;
        this.dimension = data.getDimension();
        this.firstElementPrecision = firstElementPrecision;
        this.precision = precision;
//        this.logScale = logScale;
    }

    @Override
    public double[] nextRandom() {
        double[] result = new double[dimension];
        result[0] = MathUtils.nextGaussian() / Math.sqrt(firstElementPrecision.getParameterValue(0));
        for (int i = 1; i < dimension; i++) {
            double eps = MathUtils.nextGaussian();
            result[i] = result[i-1] + eps / Math.sqrt(precision.getParameterValue(0));
        }
        return result;
    }

    @Override
    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return logPdf(v);
    }

    public double logPdf(double[] x) {
        double result = 0;

        result += NormalDistribution.logPdf(x[0], 0, 1/Math.sqrt(firstElementPrecision.getParameterValue(0)));
        for (int i = 1; i < dimension; i++) {
            result += NormalDistribution.logPdf(x[i], x[i-1], 1/Math.sqrt(precision.getParameterValue(0)));
        }

        return result;
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        return new double[0][];
    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        return logPdf(data.getParameterValues());
    }

    @Override
    public void makeDirty() {

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    public static void main(String[] args) throws Exception {
//        Parameter data = new Parameter.Default(100);
//        Parameter firstElemPrec = new Parameter.Default(0.25);
//        Parameter prec = new Parameter.Default(4.0);
//        RandomWalkGenerator gen = new RandomWalkGenerator(data, firstElemPrec, prec);
//
//        System.out.println("Data = " + Arrays.toString(data.getParameterValues()));
//
//        double[] newData = gen.nextRandom();
//
//        System.out.println("New Data = " + Arrays.toString(newData));
//
//        System.out.print("dat = c(");
//        for (int i = 0; i < newData.length; i++) {
//            data.setParameterValue(i, newData[i]);
//            System.out.printf("%.3f,", newData[i]);
//        }
//        System.out.println();
//
//        System.out.println("New Data too = " + Arrays.toString(data.getParameterValues()));
//
//        System.out.println("Log-likelihood = " + gen.getLogLikelihood());
//
//        double[] values = new double[5];
//        System.out.println("New Data = " + Arrays.toString(values));

        NewickImporter importer = new NewickImporter("((((t75_0:8.196509369,(t76_0:5.19700607,t77_0:5.117706795):2.208323646):3.611573553,(t78_0:10.09414276,(t79_0:1.555758051,t80_0:0.6394209368):8.440019033):0.7958799608):36.90916229,(t81_0:6.965230429,t82_0:6.660700811):39.60826152):84.84805599,(((((t83_0:8.461009149,(((t84_0:6.681657444,t85_0:6.329850896):0.4970932551,(t86_0:5.709749921,t87_0:5.207473638):1.032134073):0.2907208781,t88_0:5.860923226):0.6714693579):0.9286485984,(t89_0:4.35176728,(t90_0:3.383962439,t91_0:2.365432603):0.5895519359):2.156672574):0.8045815339,(((t92_0:2.347375928,(t93_0:1.806692774,t95_0:0.2737786494):0.06086363256):1.778880139,t96_0:1.745952006):0.1290466242,t97_0:1.831227907):1.231195333):2.445021787,((t98_0:0.8797524272,t99_0:0.8146248384):2.02910329,t100_0:2.379928659):2.546864865):8.217636348,(t94_0:0.1096345611,(((((t72_0:3.137047799,(t65_0:0.7428752174,(((t54_0:2.298056151,t55_0:1.386634319):0.07055842264,(t56_0:0.6211411999,t57_0:0.5106794843):0.7290213103):1.302958439,(((t32_0:0.01892218578,(((t28_0:0.8330661569,((((t25_0:0.1067751399,t1_0:9.307914765):0.2959298825,t2_0:9.425199804):0.07190576719,t3_0:8.826132669):0.2342818567,(t4_0:8.321528621,t5_0:8.0599082):0.2509419834):0.9408918394):0.7425197124,t6_0:9.694992029):2.36508096,t7_0:11.2758228):0.1994631268):8.443974525,((t8_0:10.76752205,((((t9_0:1.493044826,(t10_0:0.6698090384,t11_0:0.2989691445):0.7176777198):2.286111176,t12_0:2.654871592):1.588547853,t13_0:4.203655954):3.489759422,t14_0:7.542791391):1.88216266):8.014313382,((((t15_0:4.455846169,(t16_0:0.9348779818,t17_0:0.7092212787):2.927543424):0.8337954126,t18_0:4.371535283):2.561638646,(t19_0:2.317319668,(t20_0:0.2652308407,t21_0:0.2469742062):1.962177717):3.646397641):4.201926648,(((t22_0:3.626205298,(t23_0:1.814394049,t24_0:1.722932275):1.576199638):0.6832367053,(t26_0:0.143315467,t27_0:0.135478801):1.855008):4.430521602,(t29_0:0.277418938,t30_0:0.2751244534):2.656950561):1.277291788):4.44558088):0.4877611393):1.124425501,(((t31_0:2.199804416,t33_0:1.72439591):5.99779354,t34_0:7.697646612):0.2314130313,t35_0:7.7292073):1.578317139):0.6787721178):5.713754691):4.899357656):2.690189889,((((((t36_0:8.464440253,((t37_0:1.516370472,(t38_0:0.4518582713,t39_0:0.3776256558):1.039775741):5.008595559,(t40_0:1.662591706,t41_0:1.311689096):4.010173892):1.134735834):0.9467574058,t42_0:7.365820535):0.4308090259,((((t43_0:4.977793922,t44_0:4.1359021):0.04952307756,t45_0:4.176477429):0.004420431545,t46_0:3.828337992):0.1806606434,((t47_0:2.488665975,t48_0:2.460902709):0.2186504359,(t49_0:1.415777402,t50_0:1.358357241):0.6077191637):0.7478779126):1.81853207):7.625751786,(t51_0:1.460127624,t52_0:1.344261146):10.22672898):4.701163396,((t53_0:9.827071398,(t58_0:1.163158974,t59_0:0.8080003353):4.25790361):5.797987449,(t60_0:8.676431981,(t61_0:0.7238168631,t62_0:0.5772959095):6.677580493):1.790578221):0.5861914293):0.02660969378,t63_0:9.157389255):1.055890624):1.298920987,((t64_0:1.403696795,t66_0:0.17472493):4.761586473,(t67_0:0.8682056765,t68_0:0.7870286654):3.996084715):3.608662392):1.291247032,(t69_0:3.938693436,(t70_0:1.363501591,t71_0:1.047508982):1.675768049):5.725301901):1.748629464,(t73_0:0.4162471532,t74_0:0.3326790869):5.542586246):2.615577619):15.02704211):110.212517)");
        Tree tree = importer.importNextTree();
        Parameter N0 = new Parameter.Default(new double[] {10831.41518, 18564.20864, 4012.90013, 4844.7051, 2191.07171, 864.58823, 1954.18284, 7233.75831, 1624.29542, 1204.22108, 2424.79383, 3236.60901, 741.81022, 309.39496, 739.2772, 436.22362, 793.08349, 708.42725, 470.52578, 731.14609, 609.88306, 380.00053, 333.14087, 445.31313, 385.08516, 310.62223, 117.55287, 15.35954, 5.30181, 9.55051, 31.10367, 59.87409, 53.83162, 28.2465, 57.93228, 26.31053, 61.78336, 98.70132, 104.77869, 255.49281, 715.66912, 324.43738, 163.67341, 89.23281, 91.41586, 133.36657, 71.08849, 93.51971, 329.73119, 595.42027, 274.05989, 192.55929, 200.60026, 180.83781, 225.98434, 96.75023, 108.98978, 169.2301, 174.38823, 206.05087, 24.43085, 32.28332, 29.62218, 22.61263, 67.1912, 12.66908, 14.31542, 17.22216, 40.34572, 29.89561, 269.60226, 792.19973, 1639.95258, 495.38584, 438.89219, 266.81732, 943.54934, 1430.35456, 454.51251, 310.67618, 500.64104, 308.28057, 289.18023, 1065.39065, 633.67468, 2127.13342, 4534.01476, 1561.61864, 1958.98386, 2209.80371, 655.44965, 352.61598, 107.18735, 53.62083, 124.80736, 189.3125, 67.15406, 7.91813, 15.98184, 26.89309});
        double[] epochLengths = new double[] {1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024, 1.699024};

        DemographicModel Ne = new PiecewisePopulationModel("Ne(t)", N0, epochLengths, false, Units.Type.DAYS);

        TreeIntervals intervalList = new TreeIntervals(tree, null, null);
        CoalescentLikelihood coalescent = new CoalescentLikelihood(intervalList, Ne);

        double logLik = coalescent.getLogLikelihood();

        System.out.printf("Loglik = %f\n", logLik);

        Parameter prec = new Parameter.Default(2.0);
        Parameter logPop = new TransformedParameter(N0, new Transform.LogTransform());
        RandomWalkGenerator rwg = new RandomWalkGenerator(logPop, new Parameter.Default(0.01), prec);

        double logFieldPri = rwg.getLogLikelihood();

        System.out.printf("LogFieldPri = %f\n", logFieldPri);

//        ArrayList<Double> al = new ArrayList<Double>();
//        List<Double> lst = new Collections.synchronizedList(al);
    }

    private final Parameter data;
    private int dimension;
    private Parameter firstElementPrecision;
    private Parameter precision;
//    private final boolean logScale;
}
