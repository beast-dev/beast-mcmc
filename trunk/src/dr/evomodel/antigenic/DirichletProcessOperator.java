package dr.evomodel.antigenic;

import dr.inference.distribution.DirichletProcessLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * An integer uniform operator for allocation of items to clusters.
 *
 * @author Andrew Rambaut
 * @version $Id: UniformOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class DirichletProcessOperator extends SimpleMCMCOperator implements GibbsOperator {
    public final static String DIRICHLET_PROCESS_OPERATOR = "dirichletProcessOperator";

    private final int N;
    private final DirichletProcessLikelihood dirichletProcess;

    public DirichletProcessOperator(Parameter clusteringParameter, DirichletProcessLikelihood dirichletProcess, double weight) {
        this.clusteringParameter = clusteringParameter;
        this.N = dirichletProcess.getN();
        this.dirichletProcess = dirichletProcess;

        setWeight(weight);
    }


    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return (Parameter) clusteringParameter;
    }

    /**
     * @return the Variable this operator acts on.
     */
    public Variable getVariable() {
        return clusteringParameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        int index = MathUtils.nextInt(clusteringParameter.getDimension());

        int[] occupancy = new int[N];

        // construct cluster occupancy vector excluding the selected item and count
        // the unoccupied clusters.

        int X = N;
        for (int i = 0; i < clusteringParameter.getDimension(); i++) {
            int j = (int) clusteringParameter.getParameterValue(i);
            if (i != index) {
                occupancy[j] += 1;
                if (occupancy[j] == 1) { // first item in cluster
                    X -= 1;
                }
            }
        }

        double chi = dirichletProcess.getChiParameter().getParameterValue(0);

        double p1 = chi / ((N - 1 + chi) * X);
        double[] P = new double[N];
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            double p;
            if (occupancy[i] == 0) {
                p = p1;
            } else {
                p = occupancy[i] / (N - 1 + chi);
            }

            sum += p;
            P[i] = sum;
        }

//        if (sum != 1.0) {
//            throw new RuntimeException("doesn't sum to 1");
//        }

        double r = MathUtils.nextDouble();
        int k = 0;
        while (k < (N - 1) && r > P[k]) {
            k++;
        }

        ((Parameter) clusteringParameter).setParameterValue(index, k);

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return DIRICHLET_PROCESS_OPERATOR+"(" + clusteringParameter.getId() + "|" + dirichletProcess.getId() + ")";
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String toString() {
        return DIRICHLET_PROCESS_OPERATOR + "(" + clusteringParameter.getId() + ")";
    }

    //PRIVATE STUFF

    private Parameter clusteringParameter = null;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIRICHLET_PROCESS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            Parameter clusteringParameter = (Parameter) xo.getChild(Parameter.class);
            DirichletProcessLikelihood dirichletProcess = (DirichletProcessLikelihood) xo.getChild(DirichletProcessLikelihood.class);

            return new DirichletProcessOperator(clusteringParameter, dirichletProcess, weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
        }

        public Class getReturnType() {
            return DirichletProcessOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(DirichletProcessLikelihood.class),
                new ElementRule(Parameter.class)
        };
    };

    public int getStepCount() {
        return 1;
    }
}
