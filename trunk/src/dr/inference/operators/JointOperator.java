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

//    private ColourSamplerModel colouringModel;

    private double acceptanceFactor = 1.0; // For compatiability with ColourOperator

    private ArrayList<MCMCOperator> operator;
    private MCMCOperator innerOperator;

    public JointOperator() {

        operator = new ArrayList<MCMCOperator>();

    }

    public void addOperator(MCMCOperator operation) {
        operator.add(operation);
        if( operator.size() == 1 )
            innerOperator = operation;
    }

   public final double operate() throws OperatorFailedException {

        double logP = 0;

        for(MCMCOperator operation : operator) {
//            System.err.println("calling "+operation.getOperatorName());
//            try {
            logP += operation.operate();
//            } catch (OperatorFailedException ofe) {
//                System.err.println("catch failure in "+operation.getOperatorName()+": "+ofe.getMessage());
//                throw ofe;
//            }
//            System.err.println("operator done.");
        }

       return logP;
    }



    public double getCoercableParameter() {
        if (innerOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) innerOperator).getCoercableParameter();
        }
        throw new IllegalArgumentException();
    }

    public void setCoercableParameter(double value) {
        if (innerOperator instanceof CoercableMCMCOperator) {
            ((CoercableMCMCOperator) innerOperator).setCoercableParameter(value);
            return;
        }
        throw new IllegalArgumentException();
    }

    public double getRawParameter() {

        if (innerOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) innerOperator).getRawParameter();
        }
        throw new IllegalArgumentException();
    }

    public int getMode() {
        if (innerOperator instanceof CoercableMCMCOperator) {
            return ((CoercableMCMCOperator) innerOperator).getMode();
        }
        return CoercableMCMCOperator.COERCION_OFF;
    }

    public String getOperatorName() {
        StringBuffer sb = new StringBuffer("Joint( ");
        for(MCMCOperator operation : operator)
            sb.append(operation.getOperatorName()+" ");
        sb.append(")");
        return sb.toString();
    }

    public Element createOperatorElement(Document d) {
        throw new RuntimeException("not implemented");
    }

    public double getTargetAcceptanceProbability() {
        return innerOperator.getTargetAcceptanceProbability() * acceptanceFactor;
    }

    public double getMinimumAcceptanceLevel() {
        return innerOperator.getMinimumAcceptanceLevel() * acceptanceFactor;
    }

    public double getMaximumAcceptanceLevel() {
        return innerOperator.getMaximumAcceptanceLevel() * acceptanceFactor;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return innerOperator.getMinimumGoodAcceptanceLevel() * acceptanceFactor;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return innerOperator.getMaximumGoodAcceptanceLevel() * acceptanceFactor;
    }

    // All of this is copied and modified from SimpleMCMCOperator
    /**
     * @return the weight of this operator.
     */
    public final double getWeight() {
        return innerOperator.getWeight();
    }

    /**
     * Sets the weight of this operator.
     */
    public final void setWeight(double w) {
        innerOperator.setWeight(w);
    }

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
        return innerOperator.getAccepted();
    }

    public final void setAccepted(int accepted) {
        innerOperator.setAccepted(accepted);
    }

    public final int getRejected() {
        return innerOperator.getRejected();
    }

    public final void setRejected(int rejected) {
        innerOperator.setRejected(rejected);
    }

    public final double getMeanDeviation() {
        return innerOperator.getMeanDeviation();
    }

    public final double getSumDeviation() {
        return innerOperator.getSumDeviation();
    }

    public double getSpan(boolean reset) {
        return 0;
    }

    public final void setDumDeviation(double sumDeviation) {
        innerOperator.setDumDeviation(sumDeviation);
    }

    public String getPerformanceSuggestion() {
        return innerOperator.getPerformanceSuggestion();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return JOINT_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) {

            JointOperator operator = new JointOperator();

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

        };

    };
}

