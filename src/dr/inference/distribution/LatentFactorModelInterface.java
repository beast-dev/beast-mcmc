package dr.inference.distribution;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

public interface LatentFactorModelInterface {
    MatrixParameterInterface getFactors();

    Parameter getMissingIndicator();

    MatrixParameterInterface getLoadings();

    MatrixParameterInterface getColumnPrecision();

    MatrixParameterInterface getScaledData();

    int getFactorDimension();
}
