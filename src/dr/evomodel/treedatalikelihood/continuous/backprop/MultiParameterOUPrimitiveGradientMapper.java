package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evolution.tree.NodeRef;
import org.ejml.data.DenseMatrix64F;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class MultiParameterOUPrimitiveGradientMapper
        implements ContinuousTraitBackpropGradient.OUPrimitiveGradientMapper {

    private final List<SingleParameterOUPrimitiveGradientMapper> mappers = new ArrayList<>();
    private int totalDimension = 0;

    public void addMapper(SingleParameterOUPrimitiveGradientMapper mapper) {
        if (mapper == null) {
            throw new IllegalArgumentException("Mapper cannot be null");
        }
        mappers.add(mapper);
        totalDimension += mapper.getDimension();
    }

    public int getDimension() {
        return totalDimension;
    }

    public List<SingleParameterOUPrimitiveGradientMapper> getMappers() {
        return Collections.unmodifiableList(mappers);
    }

    @Override
    public double[] mapPrimitiveToParameters(NodeRef node,
                                             DenseMatrix64F dLdS,
                                             DenseMatrix64F dLdSigmaStat,
                                             DenseMatrix64F dLdMu,
                                             DenseMatrix64F dLdSigma) {
        double[] result = new double[totalDimension];
        int offset = 0;
        for (SingleParameterOUPrimitiveGradientMapper mapper : mappers) {
            mapper.mapPrimitiveToParameter(node, dLdS, dLdSigmaStat, dLdMu, dLdSigma, result, offset);
            offset += mapper.getDimension();
        }
        return result;
    }
}
