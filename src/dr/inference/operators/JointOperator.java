package dr.inference.operators;

import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 *
 */
public class JointOperator implements CoercableMCMCOperator {

    public static final String JOINT_OPERATOR = "jointOperator";
    public static final String WEIGHT = "weight";
    public static final String OPTIMIZE = "toOptimize";

//    private ColourSamplerModel colouringModel;

    private double acceptanceFactor = 1.0; // For compatiability with ColourOperator

    private ArrayList<MCMCOperator> operator;
    private MCMCOperator optimizedOperator;
    private int toOptimize;
    private double weight;

    public JointOperator() {
        operator = new ArrayList<MCMCOperator>();
    }

    public JointOperator(double weight, int toOptimize) {

        operator = new ArrayList<MCMCOperator>();
        this.toOptimize = toOptimize;

        setWeight(weight);

    }

    public void addOperator(MCMCOperator operation) {
        operator.add(operation);
        if( operator.size() == toOptimize )
            optimizedOperator = operation;
    }

   public final double operate() throws OperatorFailedException {

        double logP = 0;

       boolean failed = false;
       OperatorFailedException failure = null;

        for(MCMCOperator operation : operator) {
//            System.err.println("calling "+operation.getOperatorName());

                try {
                    logP += operation.operate();
                } catch (OperatorFailedException ofe) {
//                    System.err.println("catch failure in "+operation.getOperatorName()+": "+ofe.getMessage());
                    failed = true;
                    failure = ofe;
                }
//            System.err.println("operator done.");
            // todo After a failure, should not have to complete remaining operations, need to fake their operate();
        }
       if( failed )
        throw failure;

       return logP;
    }



    public double getCoercableParameter() {
        if (optimizedOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) optimizedOperator).getCoercableParameter();
        }
        throw new IllegalArgumentException();
    }

    public void setCoercableParameter(double value) {
        if (optimizedOperator instanceof CoercableMCMCOperator) {
            ((CoercableMCMCOperator) optimizedOperator).setCoercableParameter(value);
            return;
        }
        throw new IllegalArgumentException();
    }

    public double getRawParameter() {

        if (optimizedOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) optimizedOperator).getRawParameter();
        }
        throw new IllegalArgumentException();
    }

    public int getMode() {
        if (optimizedOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) optimizedOperator).getMode();
        }
        return CoercableMCMCOperator.COERCION_OFF;
    }

    public String getOperatorName() {
        StringBuffer sb = new StringBuffer("Joint(\n");
        for(MCMCOperator operation : operator)
            sb.append("\t"+operation.getOperatorName()+"\n");
        sb.append(") opt = "+optimizedOperator.getOperatorName());
        return sb.toString();
    }

    public Element createOperatorElement(Document d) {
        throw new RuntimeException("not implemented");
    }

    public double getTargetAcceptanceProbability() {
        return optimizedOperator.getTargetAcceptanceProbability() * acceptanceFactor;
    }

    public double getMinimumAcceptanceLevel() {
        return optimizedOperator.getMinimumAcceptanceLevel() * acceptanceFactor;
    }

    public double getMaximumAcceptanceLevel() {
        return optimizedOperator.getMaximumAcceptanceLevel() * acceptanceFactor;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return optimizedOperator.getMinimumGoodAcceptanceLevel() * acceptanceFactor;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return optimizedOperator.getMaximumGoodAcceptanceLevel() * acceptanceFactor;
    }

    // All of this is copied and modified from SimpleMCMCOperator

    public final void accept(double deviation) {
        for(MCMCOperator operation : operator)
            operation.accept(deviation);
    }

    public final void reject() {
        for(MCMCOperator operation : operator)
            operation.reject();
    }

    public final void reset() {
         for(MCMCOperator operation : operator)
            operation.reset();
    }

    public final int getAccepted() {
        return optimizedOperator.getAccepted();
    }

    public final void setAccepted(int accepted) {
        optimizedOperator.setAccepted(accepted);
    }

    public final int getRejected() {
        return optimizedOperator.getRejected();
    }

    public final void setRejected(int rejected) {
        optimizedOperator.setRejected(rejected);
    }

    public final double getMeanDeviation() {
        return optimizedOperator.getMeanDeviation();
    }

    public final double getSumDeviation() {
        return optimizedOperator.getSumDeviation();
    }

    public double getSpan(boolean reset) {
        return 0;
    }

    public final void setDumDeviation(double sumDeviation) {
        optimizedOperator.setDumDeviation(sumDeviation);
    }

    public String getPerformanceSuggestion() {
        return optimizedOperator.getPerformanceSuggestion();
    }

     public final double getWeight() {
        return weight;
    }

    /**
     * Sets the weight of this operator.
     */
    public final void setWeight(double w) {
        if (w > 0) {
            weight = w;
        } else throw new IllegalArgumentException("Weight must be a positive real. (called with " + w + ")");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return JOINT_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) {

            double weight;
            int toOptimize = 1;

            try {
                weight = xo.getDoubleAttribute(WEIGHT);

                if( xo.hasAttribute(OPTIMIZE) )
                    toOptimize = xo.getIntegerAttribute(OPTIMIZE);

                if( toOptimize < 0 || toOptimize > xo.getChildCount() )
                    throw new RuntimeException("Can only optimize existing operator");

            } catch (XMLParseException e) {
                throw new RuntimeException("Must provide valid 'weight' attribute");
            }

            JointOperator operator = new JointOperator(weight,toOptimize);

            for(int i=0; i<xo.getChildCount(); i++) {
                operator.addOperator((MCMCOperator)xo.getChild(i));

            }

           return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an arbitrary list of operators; only the first is optimizable";
        }

        public Class getReturnType() {
            return JointOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(MCMCOperator.class,1, Integer.MAX_VALUE),
                AttributeRule.newDoubleArrayRule(WEIGHT),
                AttributeRule.newIntegerArrayRule(OPTIMIZE,true)
        };

    };
}

