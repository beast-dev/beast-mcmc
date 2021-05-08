package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.function.DoubleBinaryOperator;

public class RandomLocalClockRateWrapper implements DifferentiableBranchRates {

    public static final String RLC_RATES_WRAPPER = "rlcRatesWrapper";
//    public static final String PARAPHYLY_LIST = "paraphylyList";
//    public static final String WEIGHTING = "weighting";

    //    private List<Taxa> paraphylySet;
    private RandomLocalClockModel rlcModel;
//    private Tree tree;
//    private double totalTime;
//    private List<NodeRef> MRCANodeList;
//    private int dim;
//    private BranchWeighting branchWeighting;


    public RandomLocalClockRateWrapper(RandomLocalClockModel rlcModel) {

        this.rlcModel = rlcModel;
//        super(name);
//        this.branchRateModel = branchRateModel;
//        this.paraphylySet = paraphylySet;
//        this.tree = branchRateModel.getTree();
//        this.branchWeighting = branchWeighting;
//        this.dim = dim;
//        this.MRCANodeList = new ArrayList<NodeRef>(dim);
//        for (int i = 0; i < dim; i++) {
//            MRCANodeList.add(null);
//        }
    }

    @Override
    public double getUntransformedBranchRate(Tree tree, NodeRef node) {
//        return getBranchRate(tree, node);
        return rlcModel.getUnscaledBranchRate(tree, node);
    }


    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getRateParameter() {
        return null;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return 0;
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        return null;
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");

    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");

    }

    @Override
    public double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        return 0;
    }

    @Override
    public void forEachOverRates(NodeRateMap map) {

    }

    @Override
    public double getPriorRateAsIncrement(Tree tree) {
        return 0;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        return 0;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return RLC_RATES_WRAPPER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//            final String name = xo.getAttribute(Statistic.NAME, xo.getId());

//            List<String> names = new ArrayList<>();

            RandomLocalClockModel rlcModel = (RandomLocalClockModel) xo.getChild(RandomLocalClockModel.class);

//            List<Taxa> paraphylySet = new ArrayList<>();
//
//            Taxa paraphyly;
//
//            if (xo.hasChildNamed(PARAPHYLY_LIST)) {
//                XMLObject cxo = xo.getChild(PARAPHYLY_LIST);
//                for (int i = 0; i < cxo.getChildCount(); i++) {
//                    paraphyly = (Taxa) cxo.getChild(i);
//                    paraphylySet.add(paraphyly);
//                }
//            }
//
//            int dim = paraphylySet.size();
//
//            BranchWeighting branchWeighting = parseWeighting(xo);

            RandomLocalClockRateWrapper rlcRateWrapper = new RandomLocalClockRateWrapper(rlcModel);
            return rlcRateWrapper;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RandomLocalClockModel.class),
        };

        public String getParserDescription() {
            return "DifferentiableBranchRates object that wraps around randomLocalClock rates (product of rates and indicators)";
        }

        public Class getReturnType() {
            return RandomLocalClockRateWrapper.class;
        }
    };

}
