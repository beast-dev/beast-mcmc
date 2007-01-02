package dr.inference.distribution;

import dr.inference.loggers.LogColumn;
import dr.inference.model.*;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Dec 29, 2006
 * Time: 11:01:02 AM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralizedLinearModel extends AbstractModel implements Likelihood {


    public static final String GLM_LIKELIHOOD = "glmLikelihood";

    public static final String DEPENDENT_VARIABLES = "dependentVariables";
    public static final String INDEPENDENT_VARIABLES = "independentVariables";
    public static final String BASIS_MATRIX = "basis";

    private Parameter dependentParam;
    private Parameter independentParam;
    private DesignMatrix designMatrix;

    public GeneralizedLinearModel(Parameter dependentParam, Parameter independentParam,
                                  DesignMatrix designMatrix) {
        super(GLM_LIKELIHOOD);
        this.dependentParam = dependentParam;
        this.independentParam = independentParam;
        this.designMatrix = designMatrix;
        addParameter(independentParam);
        addParameter(dependentParam);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {
    }

    protected void storeState() {
        // No internal states to save

    }

    protected void restoreState() {
        // No internal states to restore
    }

    protected void acceptState() {
        // Nothing to do
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return calculateLogLikelihood();
//        return 0;
    }


    @Override
    public String toString() {
        return super.toString() + ": " + getLogLikelihood();
    }

    private double calculateLogLikelihood() {
        // logLikelihood calculation for logistic regression
        double logLikelihood = 0;

        final int K = independentParam.getDimension();
        final int N = dependentParam.getDimension();

        for (int i = 0; i < N; i++) {
            // for each "pseudo"-datum
            double xBeta = 0;
            for (int k = 0; k < K; k++) {
                xBeta += designMatrix.getParameterValue(i, k) * independentParam.getParameterValue(k);
            }

            logLikelihood += dependentParam.getParameterValue(i) * xBeta
                    - Math.log(1.0 + Math.exp(xBeta));

        }
        return logLikelihood;
    }

    public void makeDirty() {
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }


    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GLM_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(DEPENDENT_VARIABLES);
            Parameter dependentParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(INDEPENDENT_VARIABLES);
            Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(BASIS_MATRIX);
            DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);

            if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
                    (independentParam.getDimension() != designMatrix.getColumnDimension()))
                throw new XMLParseException(
                        "dim(" + DEPENDENT_VARIABLES + ") != dim(" + BASIS_MATRIX + " %*% " + INDEPENDENT_VARIABLES + ")"
                );

            return new GeneralizedLinearModel(dependentParam, independentParam, designMatrix);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(DEPENDENT_VARIABLES,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(INDEPENDENT_VARIABLES,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(BASIS_MATRIX,
                        new XMLSyntaxRule[]{new ElementRule(DesignMatrix.class)})
//                new ElementRule(COUNTS,
//                        new XMLSyntaxRule[]{AttributeRule.newIntegerArrayRule("values", false),})
        };

        public String getParserDescription() {
            return "Calculates the generalized linear model likelihood of the dependent parameters given the indepenent parameters and design matrix.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };
}
