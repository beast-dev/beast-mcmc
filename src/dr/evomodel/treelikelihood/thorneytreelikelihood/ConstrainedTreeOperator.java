package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evomodel.operators.AbstractAdaptableTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import dr.math.Poisson;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

public class ConstrainedTreeOperator extends AbstractAdaptableTreeOperator {
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String MEAN_COUNT = "meanCount";
    public static final String FIXED_COUNT = "fixedCount";

    public static ConstrainedTreeOperator parse(ConstrainedTreeModel tree, double weight, ConstrainableTreeOperator op, XMLObject xo) throws XMLParseException {
        AdaptationMode mode = AdaptationMode.parseMode(xo);

        int count = xo.getAttribute(FIXED_COUNT,1);
        double nonShiftedMean = xo.getAttribute(MEAN_COUNT, 1.0);
        double targetAcceptanceProbability = xo.getAttribute("targetAcceptance", 0.234D);
        if(!xo.hasAttribute(MEAN_COUNT)){
            mode=AdaptationMode.ADAPTATION_OFF;
        }
        if(xo.hasAttribute(MEAN_COUNT) && xo.hasAttribute(FIXED_COUNT)){
            throw new XMLParseException("Constrained operator with id " + op.getOperatorName()+":  meanCount and fixedCount are mutually exclusive");
        }
        return new ConstrainedTreeOperator(tree,weight,op,nonShiftedMean,count,mode,targetAcceptanceProbability);
    }

    public ConstrainedTreeOperator(ConstrainedTreeModel tree, double weight, ConstrainableTreeOperator operator,double nonshiftedMean,int count,AdaptationMode mode, double targetAcceptanceProbability) {
        super(mode, targetAcceptanceProbability);
        setWeight(weight);
        constrainedTreeModel = tree;
        this.operator = operator;
        if (tree.getInternalNodeCount() == tree.getSubtreeCount()) {
            throw new IllegalArgumentException(getOperatorName() + " is designed to resolve polytomies; however, the "+
                    "constrained tree is fully resolved. Please remove this operator or provide an unresolved "+
                    "constraints tree.");
        }
        subtreeSizes = new double [tree.getSubtreeCount()];
        maxOperations=0;
        for (int i = 0; i < tree.getSubtreeCount(); i++) {
            subtreeSizes[i] =  (double) tree.getSubtree(i).getInternalNodeCount()-1; // don't choose subtrees with only 1 internal node. There's no topology to sample.
            if(subtreeSizes[i]>0){
                maxOperations++;
            }
        }

        this.nonshiftedMean = nonshiftedMean;
        this.count=count;
    }

    @Override
    public String getOperatorName() {
        return "Constrained  " + operator.getOperatorName();
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    @Override
    public double doOperation() {

        int operations = getOperationCount();
        double[] selectionSizes = new double[subtreeSizes.length];
        System.arraycopy(subtreeSizes,0,selectionSizes,0,subtreeSizes.length);
        int[] selected = new int[operations];

        for (int i = 0; i < operations; i++) {
            int selection =  MathUtils.randomChoicePDF(selectionSizes);
            selected[i] =  selection;
            selectionSizes[selection] = 0.0;
        }

            double logP = 0;
            for (int i =0;i<operations;i++) {
                TreeModel subtree = constrainedTreeModel.getSubtree(selected[i]);
                logP += operator.doOperation(subtree);
            }
            return logP;
        }


        /**
         * Sets the adaptable parameter value.
         *
         * @param value the value to set the adaptable parameter to
         */
        @Override
        protected void setAdaptableParameterValue(double value) {
            nonshiftedMean = Math.exp(value);
        }

        /**
         * Gets the adaptable parameter value.
         *
         * @returns the value
         */
        @Override
        protected double getAdaptableParameterValue() {
            return Math.log(nonshiftedMean);
        }

        /**
         * @return the underlying tuning parameter value
         */
        @Override
        public double getRawParameter() {
            return nonshiftedMean;
        }

        @Override
        public String getAdaptableParameterName() {
            return "unshiftedMean";
        }

        private int getOperationCount(){
            if(this.getMode()== AdaptationMode.ADAPTATION_OFF){
                return this.count;
            }else{
                int operations = Poisson.nextPoisson(this.nonshiftedMean) +1;
                return Math.min(operations, maxOperations);
            }
    }

    private final ConstrainedTreeModel constrainedTreeModel;
    private final ConstrainableTreeOperator operator;
    private final double[] subtreeSizes;
    private double nonshiftedMean;
    private final int count;
    private  int maxOperations;


}
