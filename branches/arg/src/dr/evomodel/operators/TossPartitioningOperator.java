package dr.evomodel.operators;

import dr.evomodel.tree.ARGModel;
import dr.inference.model.Parameter;
import dr.inference.model.VariableSizeCompoundParameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TossPartitioningOperator extends SimpleMCMCOperator {

    private VariableSizeCompoundParameter partitioningParameters;
    private ARGModel argModel;
    //private int numberPartitions;

    private final static String OPERATOR_NAME = "tossPartitioningOperator";

    public TossPartitioningOperator(ARGModel arg) {
        argModel = arg;
        partitioningParameters = arg.getPartitioningParameters();
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return partitioningParameters;
    }

    @Override
    public final double doOperation() throws OperatorFailedException {
        int numberParameters = partitioningParameters.getNumParameters();
        if (numberParameters == 0)
            throw new OperatorFailedException("Must be atleast one reassortment.");
        // Pick one
        /* else     {
if (true)

throw new OperatorFailedException("out");         }*/
        // TODO update more than one at a time
        int whichParameter = MathUtils.nextInt(numberParameters);
        Parameter partitioning = partitioningParameters.getParameter(whichParameter);
        int numberPartitions = partitioning.getDimension();
	    //System.err.println("np = "+numberPartitions);
	    //System.exit(-1);
        if (numberPartitions == 2) { // Just swap partitioning
            if (partitioning.getParameterValue(0) == 0) {
                partitioning.setParameterValueQuietly(0, 1);
                partitioning.setParameterValueQuietly(1, 0);
            } else {
                partitioning.setParameterValueQuietly(0, 0);
                partitioning.setParameterValueQuietly(1, 1);
            }
        } else { // There are more than just two possible partitioning
            // generate a uniform random draw from all possible partitionings
            int[] permutation = MathUtils.permuted(numberPartitions);
            int cut = MathUtils.nextInt(numberPartitions - 1);
            for (int i = 0; i < numberPartitions; i++) {
                if (i <= cut)
                    partitioning.setParameterValueQuietly(permutation[i], 0);
                else
                    partitioning.setParameterValueQuietly(permutation[i], 1);
            }

        }
        argModel.fireModelChanged(new PartitionChangedEvent(partitioningParameters));
        return 0;
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

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//			int weight = xo.getIntegerAttribute(WEIGHT);
//			Parameter parameter = (Parameter)xo.getChild(Parameter.class);	
            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);
            return new TossPartitioningOperator(arg);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that picks a new partitioning uniformly at random.";
        }

        public Class getReturnType() {
            return TossPartitioningOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//			AttributeRule.newIntegerRule(WEIGHT),
                new ElementRule(ARGModel.class)
        };
    };

    public String toString() {
        return "tossPartitioningOperator(" + partitioningParameters.getParameterName() + ")";
    }
	

}
