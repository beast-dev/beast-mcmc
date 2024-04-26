package dr.inference.operators.factorAnalysis;

import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.AbstractAdaptableOperator;
import dr.inference.operators.GeneralOperator;
import dr.inference.operators.PathDependent;
import dr.xml.*;

public class LoadingsRotationOperator extends AbstractAdaptableOperator implements GeneralOperator, PathDependent {

    private final AbstractAdaptableOperator baseOperator;
    private final MatrixParameterInterface parameter;

    public LoadingsRotationOperator(AbstractAdaptableOperator baseOperator,
                                    MatrixParameterInterface parameter
    ) {
        super(baseOperator.mode, baseOperator.getTargetAcceptanceProbability());
        this.baseOperator = baseOperator;
        this.parameter = parameter;
    }

    @Override
    public String getOperatorName() {
        return PARSER_NAME + ":" + baseOperator.getOperatorName();
    }

    @Override
    public double doOperation(Likelihood joint) {

        syncBaseOperator();

        double hastingsRatio = baseOperator.doOperation(joint);

        double oldLikelihood = joint.getLogLikelihood();

        reflect();

        double newLikelihood = joint.getLogLikelihood();

        if (oldLikelihood != newLikelihood) {
            throw new RuntimeException("Fix this");
        }

        return hastingsRatio;

    }

    private void reflect() {
        boolean changed = false;
        int dim = Math.min(parameter.getRowDimension(), parameter.getColumnDimension());
        for (int i = 0; i < dim; i++) {

            if (parameter.getParameterValue(i, i) < 0) {
                changed = true;
                for (int j = i; j < parameter.getRowDimension(); j++) {
                    parameter.setParameterValueQuietly(j, i, -parameter.getParameterValue(j, i));
                }
            }
        }

        if (changed) {
            parameter.fireParameterChangedEvent();
        }
    }

    @Override
    public double doOperation() {

        syncBaseOperator();

        double hastingsRatio = baseOperator.doOperation();
        reflect();

        return hastingsRatio;
    }

    private void syncBaseOperator() {
        baseOperator.setAcceptCount(getAcceptCount());
        baseOperator.setRejectCount(getRejectCount());
        baseOperator.setSumDeviation(getSumDeviation());
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        baseOperator.setAdaptableParameter(value);
    }

    @Override
    protected double getAdaptableParameterValue() {
        return baseOperator.getAdaptableParameter();
    }

    @Override
    public double getRawParameter() {
        return baseOperator.getRawParameter();
    }

    @Override
    public String getAdaptableParameterName() {
        return baseOperator.getAdaptableParameterName();
    }


    public static final String PARSER_NAME = "loadingsRotationOperator";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) {

            AbstractAdaptableOperator baseOperator = (AbstractAdaptableOperator)
                    xo.getChild(AbstractAdaptableOperator.class);
            MatrixParameterInterface parameter = (MatrixParameterInterface)
                    xo.getChild(MatrixParameterInterface.class);

            return new LoadingsRotationOperator(baseOperator, parameter);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(AbstractAdaptableOperator.class),
                    new ElementRule(MatrixParameterInterface.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return LoadingsRotationOperator.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };
}
