package dr.inference.model;


import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Gabriel Hassler
 * @author Marc Suchard
 */


public class ScaledParameter extends ProductParameter {

    private final Parameter scaleParam;
    private final Parameter vecParam;

    public ScaledParameter(Parameter scaleParam, Parameter vecParam) {
        super(new ArrayList<>(Arrays.asList(scaleParam, vecParam)));
        this.scaleParam = scaleParam;
        this.vecParam = vecParam;

    }

    @Override
    public double getParameterValue(int dim) {
        return scaleParam.getParameterValue(0) * vecParam.getParameterValue(dim);
    }

    @Override
    public int getDimension() {
        return vecParam.getDimension();
    }

}
