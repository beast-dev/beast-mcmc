package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * Operator that selects a microsatellite model from a group provided by the user.
 */
public class MicrosatelliteModelSelectOperator extends SimpleMCMCOperator{
    private Parameter parameter;
    private Parameter[] indicators;
    public MicrosatelliteModelSelectOperator(Parameter parameter, Parameter[] indicators, double weight){
        this.parameter = parameter;
        this.indicators = indicators;
        setWeight(weight);
    }

    public String getOperatorName(){
        return "msatModelSelectOperator("+parameter.getParameterName()+")";
    }

    public final String getPerformanceSuggestion() {
        return "no suggestions available";
    }

    public double doOperation(){
        int index = MathUtils.nextInt(indicators.length);
        //System.out.println(index);
        Parameter newModel = indicators[index];
        for(int i = 0; i < parameter.getDimension() -1 ; i++){
            parameter.setParameterValueQuietly(i,newModel.getParameterValue(i));
        }
        parameter.setParameterValueNotifyChangedAll(
                parameter.getDimension()-1,
                newModel.getParameterValue(parameter.getDimension()-1)
        );
        //System.out.println(parameter+"; "+indicators[index]);
        return 0.0;
    }



}
