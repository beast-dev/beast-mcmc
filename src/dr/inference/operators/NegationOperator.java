package dr.inference.operators;


import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * Created by maxryandolinskytolkoff on 8/17/16.
 */
public class NegationOperator extends SimpleMCMCOperator {
    Parameter data = null;

    public NegationOperator(Parameter data, Double weight){
        setWeight(weight);
        this.data = data;

    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "NegationOperator";
    }

    @Override
    public double doOperation() {
        int number = MathUtils.nextInt(data.getDimension());
        double setTo = -data.getParameterValue(number);
        data.setParameterValue(number, setTo);

        return 0;
    }
}
