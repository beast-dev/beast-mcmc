package dr.evomodel.coalescent.operators;

import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**

 */
public class SampleNonActiveGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private final ParametricDistributionModel distribution;
    private final Parameter data;
    private final Parameter indicators;

    public SampleNonActiveGibbsOperator(ParametricDistributionModel distribution,
                                        Parameter data, Parameter indicators, double weight) {
        this.distribution = distribution;
        this.data = data;
        this.indicators = indicators;
        setWeight(weight);
    }


    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return "SampleNonActive(" + indicators.getId() + ")";
    }

    public double doOperation() throws OperatorFailedException {
        final int idim = indicators.getDimension();

        final int offset = (data.getDimension() - 1) == idim ? 1 : 0;
        assert offset == 1 || data.getDimension() == idim : "" + idim + " (?+1) != " + data.getDimension();

        // available locations for direct sampling
        int[] loc = new int[idim];
        int nLoc = 0;

        for (int i = 0; i < idim; ++i) {
            final double value = indicators.getStatisticValue(i);
            if (value == 0) {
                loc[nLoc] = i + offset;
                ++nLoc;
            }
        }

        if (nLoc > 0) {
            final int index = loc[MathUtils.nextInt(nLoc)];
            try {
                final double val = distribution.quantile(MathUtils.nextDouble());
                data.setParameterValue(index, val);
            } catch (Exception e) {
                // some distributions fail on extreme values - currently gamma
               throw new OperatorFailedException(e.getMessage());
            }
        } else {
            throw new OperatorFailedException("no non-active indicators");
        }
        return 0;
    }

    public int getStepCount() {
        return 0;
    }
}
