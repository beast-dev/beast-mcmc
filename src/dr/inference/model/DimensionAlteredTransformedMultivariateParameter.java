package dr.inference.model;

import dr.util.Transform;

public class DimensionAlteredTransformedMultivariateParameter extends TransformedMultivariateParameter {

    public DimensionAlteredTransformedMultivariateParameter(Parameter parameter, Transform.MultivariableTransform transform) {
        this(parameter, transform, false);
    }

    public DimensionAlteredTransformedMultivariateParameter(Parameter parameter, Transform.MultivariableTransform transform, boolean inverse) {
        super(parameter, transform, inverse);
    }

    public int getDimension() {
        return ((Transform.MultivariableTransform) transform).getDimension();
    }

}
