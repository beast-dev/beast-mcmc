package dr.inference.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.model.Parameter;
import dr.util.Transform;

public class NewTransformedParameterOperator extends AbstractAdaptableOperator {
    private final boolean isAdaptable;
    private final SimpleMCMCOperator subOperator;
    private final Parameter parameter;
    private final Transform transform;
    private final boolean inverse;
    private final boolean checkValid;
    private final BoundedSpace generalBounds;
    public static final String NEW_TRANSFORMED_OPERATOR = "newTransformedParameterOperator";

    public NewTransformedParameterOperator(SimpleMCMCOperator operator,
                                           Transform transform,
                                           boolean inverse,
                                           BoundedSpace generalBounds) {

        this.subOperator = operator;
        this.transform = transform;
        this.inverse = inverse;
        setWeight(operator.getWeight());
        this.isAdaptable = operator instanceof AbstractAdaptableOperator;
        this.parameter = operator.getParameter();
        this.generalBounds = generalBounds;
        this.checkValid = generalBounds != null;
    }

    public NewTransformedParameterOperator(SimpleMCMCOperator operator,
                                           Transform transform,
                                           BoundedSpace generalBounds) {
        this(operator, transform, false, generalBounds);
    }


    @Override
    protected void setAdaptableParameterValue(double value) {
        if (isAdaptable) {
            ((AbstractAdaptableOperator) subOperator).setAdaptableParameterValue(value);
        }
    }

    @Override
    protected double getAdaptableParameterValue() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterValue();
        }
        return 0;
    }

    @Override
    public double getRawParameter() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getRawParameter();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getAdaptableParameterName() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterName();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getOperatorName() {
        return NEW_TRANSFORMED_OPERATOR + "." + subOperator.getOperatorName();
    }

    @Override
    public double doOperation() {
        double[] oldValues = parameter.getParameterValues();
        double ratio = subOperator.doOperation();
        double[] newValues = parameter.getParameterValues();


        if (checkValid) { // GH: below is sloppy, but best I could do without refactoring how Parameter handles bounds
            if (generalBounds == null && !parameter.isWithinBounds()) {
                return Double.NEGATIVE_INFINITY;
            } else if (!generalBounds.isWithinBounds(parameter.getParameterValues())) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Compute Jacobians

        if (inverse) {
            ratio += -transform.getLogJacobian(transform.transform(oldValues, 0, oldValues.length), 0, oldValues.length)
                    + transform.getLogJacobian(transform.transform(newValues, 0, newValues.length), 0, oldValues.length);
        } else {
            ratio += transform.getLogJacobian(oldValues, 0, oldValues.length)
                    - transform.getLogJacobian(newValues, 0, newValues.length);
        }

        return ratio;
    }

    @Override
    public Parameter getParameter() {
        return subOperator.getParameter();
    }
}
