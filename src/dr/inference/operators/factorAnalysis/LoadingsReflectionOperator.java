package dr.inference.operators.factorAnalysis;

import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.GeneralOperator;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

public class LoadingsReflectionOperator extends SimpleMCMCOperator {

    private final SimpleMCMCOperator operator;
    private final MatrixParameterInterface loadings;

    LoadingsReflectionOperator(SimpleMCMCOperator operator, MatrixParameterInterface loadings) {

        this.operator = operator;
        this.loadings = loadings;
        double weight = operator.getWeight();
        this.setWeight(weight);

    }

    @Override
    public String getOperatorName() {
        return null;
    }

    @Override
    public double doOperation() {
        double hastingsRatio = operator.doOperation();
        reflectLoadings();
        return hastingsRatio;
    }

    @Override
    public double doOperation(Likelihood likelihood) {
        double hastingsRatio = operator.doOperation(likelihood);
        reflectLoadings();
        return hastingsRatio;
    }

    private void reflectLoadings() {
        boolean changed = false;
        int dim = Math.min(loadings.getRowDimension(), loadings.getColumnDimension());
        for (int i = 0; i < dim; i++) {

            if (loadings.getParameterValue(i, i) < 0) {
                changed = true;
                for (int j = i; j < loadings.getColumnDimension(); j++) {
                    loadings.setParameterValueQuietly(i, j, -loadings.getParameterValue(i, j));
                }
            }
        }

        if (changed) {
            loadings.fireParameterChangedEvent();
        }
    }

    public static final String LOADINGS_REFLECTION_OP = "loadingsReflectionOperator";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            SimpleMCMCOperator operator = (SimpleMCMCOperator) xo.getChild(SimpleMCMCOperator.class);
            MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

            if (operator instanceof GeneralOperator && operator instanceof GibbsOperator) {
                System.err.println("Not yet implemented");
            }

            if (operator instanceof GeneralOperator) {
                return new GeneralLoadingsReflectionOperator(operator, loadings);
            } else if (operator instanceof GibbsOperator) {
                return new GibbsLoadingsReflectionOperator(operator, loadings);
            } else {
                return new LoadingsReflectionOperator(operator, loadings);
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MCMCOperator.class),
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Operator that enforces reflection constraint on child operators";
        }

        @Override
        public Class getReturnType() {
            return LoadingsReflectionOperator.class;
        }

        @Override
        public String getParserName() {
            return LOADINGS_REFLECTION_OP;
        }
    };
}
