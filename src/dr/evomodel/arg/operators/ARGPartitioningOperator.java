package dr.evomodel.arg.operators;

import dr.evomodel.arg.ARGModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.logging.Logger;

public class ARGPartitioningOperator extends SimpleMCMCOperator {

    private CompoundParameter partitioningParameters;
    private ARGModel arg;

    public final static String OPERATOR_NAME = "argPartitionOperator";
    public static final String TOSS_SIZE = "tossSize";

    private boolean isRecombination;
    private int tossSize;

    public ARGPartitioningOperator(ARGModel arg, int tossSize, int weight) {
        super.setWeight(weight);

        this.arg = arg;
        this.partitioningParameters = arg.getPartitioningParameters();
        this.tossSize = tossSize;
        this.isRecombination = arg.isRecombinationPartitionType();
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return partitioningParameters;
    }


    public final double doOperation() throws OperatorFailedException {
        double logq = 0;

        int len = partitioningParameters.getNumberOfParameters();

        if (len == 0) {
            return 0;
        }

        if (isRecombination) {
            logq = doRecombination(partitioningParameters.getParameter(MathUtils.nextInt(len)));
        } else {
            logq = doReassortment(partitioningParameters.getParameter(MathUtils.nextInt(len)));
        }

        arg.fireModelChanged(new PartitionChangedEvent(partitioningParameters));
        return logq;
    }


    private double doRecombination(Parameter partition) throws OperatorFailedException {

        assert checkValidRecombinationPartition(partition);

        int currentBreakLocation = 0;
        for (int i = 0, n = arg.getNumberOfPartitions(); i < n; i++) {
            if (partition.getParameterValue(i) == 1) {
                currentBreakLocation = i;
                break;
            }
        }

        assert currentBreakLocation > 0;

        if (MathUtils.nextBoolean()) {
            //Move break right 1
            partition.setParameterValueQuietly(currentBreakLocation, 0.0);
        } else {
            partition.setParameterValueQuietly(currentBreakLocation - 1, 1.0);
        }

        if (!checkValidRecombinationPartition(partition)) {
            throw new OperatorFailedException("");
        }


        return 0;
    }

    public static boolean checkValidRecombinationPartition(Parameter partition) {
        int l = partition.getDimension();
        if ((partition.getParameterValue(0) == 0 && partition.getParameterValue(l - 1) == 1))
            return true;

        return false;
    }


    private double doReassortment(Parameter partition) throws OperatorFailedException {

        assert checkValidReassortmentPartition(partition);

        ArrayList<Integer> list = new ArrayList<Integer>(tossSize);

        while (list.size() < tossSize) {
            int a = MathUtils.nextInt(arg.getNumberOfPartitions() - 1) + 1;
            if (!list.contains(a)) {
                list.add(a);
            }
        }


        for (int a : list) {
            if (partition.getParameterValue(a) == 0) {
                partition.setParameterValueQuietly(a, 1);
            } else {
                partition.setParameterValueQuietly(a, 0);
            }
        }
        
        
        if (!checkValidReassortmentPartition(partition)) {
            throw new OperatorFailedException("");
        }

        return 0;
    }

    public static boolean checkValidReassortmentPartition(Parameter partition) {
        if (partition.getParameterValue(0) != 0)
            return false;

        double[] a = partition.getParameterValues();

        double sum = 0;

        for (double b : a)
            sum += b;

        if (sum == 0 || sum == a.length)
            return false;

        return true;

    }

    @Override
    public String getOperatorName() {
        return OPERATOR_NAME;
    }

    public String getPerformanceSuggestion() {
        // TODO Auto-generated method stub
        return null;
    }

    public class PartitionChangedEvent {
        Parameter partitioning;

        public PartitionChangedEvent(Parameter partitioning) {
            this.partitioning = partitioning;
        }

        public Parameter getParameter() {
            return partitioning;
        }
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public String[] getParserNames() {
            return new String[]{
                    OPERATOR_NAME,
                    "tossPartitioningOperator",
            };
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int weight = xo.getIntegerAttribute(WEIGHT);

            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            int tossSize = 1;
            if (xo.hasAttribute(TOSS_SIZE)) {
                tossSize = xo.getIntegerAttribute(TOSS_SIZE);

                if (tossSize <= 0 || tossSize >= arg.getNumberOfPartitions()) {
                    throw new XMLParseException("Toss size is incorrect");
                }
            }

            Logger.getLogger("dr.evomodel").info("Creating ARGPartitionOperator with " + TOSS_SIZE + " of " + tossSize);


            return new ARGPartitioningOperator(arg, tossSize, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that picks a new partitioning uniformly at random.";
        }

        public Class getReturnType() {
            return ARGPartitioningOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(WEIGHT),
                new ElementRule(ARGModel.class)
        };
    };

    public String toString() {
        return "tossPartitioningOperator(" + partitioningParameters.getParameterName() + ")";
    }

    private boolean singleFlip;


}
