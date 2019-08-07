package dr.inference.operators.factorAnalysis;

import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.xml.*;

public class LoadingsRotationHMCOperator extends HamiltonianMonteCarloOperator {


    private final MatrixParameterInterface parameter;

    public LoadingsRotationHMCOperator(HamiltonianMonteCarloOperator hmcOp, MatrixParameterInterface parameter) {

        super(hmcOp);
        this.parameter = parameter;
    }

    @Override
    public double doOperation(Likelihood joint) {
        double hastingsRatio = super.doOperation(joint);

        boolean changed = false;
        int dim = Math.min(parameter.getRowDimension(), parameter.getColumnDimension());
        for (int i = 0; i < dim; i++) {

            if (parameter.getParameterValue(i, i) < 0) {
                changed = true;
                for (int j = i; j < parameter.getColumnDimension(); j++) {
                    parameter.setParameterValueQuietly(i, j, -parameter.getParameterValue(i, j));
                }
            }
        }

        if (changed) {
            parameter.fireParameterChangedEvent();
        }

        return hastingsRatio;

    }

    public static final String PARSER_NAME = "loadingsRotationHMCOperator";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            HamiltonianMonteCarloOperator hmcOp = (HamiltonianMonteCarloOperator)
                    xo.getChild(HamiltonianMonteCarloOperator.class);
            MatrixParameterInterface parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            return new LoadingsRotationHMCOperator(hmcOp, parameter);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(HamiltonianMonteCarloOperator.class),
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return LoadingsRotationHMCOperator.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };


}
