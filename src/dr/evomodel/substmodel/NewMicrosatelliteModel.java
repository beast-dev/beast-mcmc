package dr.evomodel.substmodel;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.evolution.datatype.Microsatellite;
import dr.math.ModifiedBesselFirstKind;

/**
 * @author Chieh-Hsi Wu
 * Implementation of models by Watkins (2007)
 */
public class NewMicrosatelliteModel extends AbstractSubstitutionModel {

    Parameter biasConst;
    public NewMicrosatelliteModel(Microsatellite msat, FrequencyModel rootFreqModel){
        super("NewMicrosatelliteModel", msat, rootFreqModel);
        biasConst = new Parameter.Default(0.5);
    }

    public void getTransitionProbabilities(double distance, double[] matrix){
        int k = 0;
        double[] rowSums = new double[stateCount];
        double bCVal = biasConst.getParameterValue(0);
        for(int i = 0; i < stateCount; i ++){
            for(int j = 0; j < stateCount; j++){
                int n = Math.abs(i - j);
                //matrix[k] = Math.exp(-distance)* Math.pow(bCVal/(1-bCVal),n/2)*ModifiedBesselFirstKind.bessi(2*Math.sqrt(bCVal*(1-bCVal))*distance,n);
                matrix[k] = Math.exp(-distance)*ModifiedBesselFirstKind.bessi(distance,n);

                rowSums[i] += matrix[k];
                k++;
            }

        }
        k = 0;
        for(int i = 0; i < stateCount; i ++){
            for(int j = 0; j < stateCount; j++){
                matrix[k] =  matrix[k]/rowSums[i];
                k++;
            }          
        }
    }

    protected void ratesChanged() {};
    protected void setupRelativeRates(){};
    protected void frequenciesChanged() {};

    public static void main(String[] args){
        Microsatellite msat = new Microsatellite(1,30);
        NewMicrosatelliteModel nmsatModel = new NewMicrosatelliteModel(msat, null);
        double[] probs = new double[msat.getStateCount()*msat.getStateCount()];
        nmsatModel.getTransitionProbabilities(1.2,probs);
        int k =0;
        for(int i = 0; i < msat.getStateCount(); i++){
            for(int j = 0; j < msat.getStateCount(); j++){
                System.out.print(probs[k++]+" ");
            }
            System.out.println();
        }
    }

}
