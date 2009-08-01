package dr.evomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.inference.model.Model;


/**
 * @author Chieh-Hsi Wu
 * An abstract class for microsatellite models
 */
public abstract class MicrosatelliteModel extends ComplexSubstitutionModel{
    protected MicrosatelliteModel subModel;

    /**
     * Constructor
     * @param name              Model name
     * @param microsatellite    Microsatellite data type
     * @param rootFreqModel     Frequency model
     * @param parameter         Infinitesimal rates
     */
    public MicrosatelliteModel(String name, Microsatellite microsatellite, FrequencyModel rootFreqModel, Parameter parameter) {
        super(name, microsatellite, rootFreqModel, parameter);
        if(parameter == null){
            double[] q = new double[stateCount*(stateCount-1)];
            infinitesimalRates = new Parameter.Default(q);
        }
    }
    


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        updateMatrix = true;
    }

    //store the infinitesimal rates in the vector to a matrix called amat
    public void storeIntoAmat(double[] rates){
        int i, j, k = 0;
        // Set the instantaneous rate matrix
        for (i = 0; i < stateCount; i++) {
            for (j = 0; j < stateCount; j++) {
                if (i != j)
                    amat[i][j] = rates[k++];
            }
        }
    }

    //setup empirical frequencies
    public void setupEmpiricalStationaryFrequencies(){
        setToEqualFrequencies();
        double[] transPrs = new double[stateCount*stateCount];
        getTransitionProbabilities(Double.MAX_VALUE, transPrs);
        double[] empFreq = new double[stateCount];
        System.arraycopy(transPrs, 0, empFreq, 0, stateCount);
        this.freqModel = new FrequencyModel(dataType, empFreq);
    }

    private void setToEqualFrequencies(){
        double[] freqs = new double[stateCount];
            for(int i = 0; i < freqs.length; i++){
                freqs[i] = 1.0/stateCount;
            }
        this.freqModel = new FrequencyModel(dataType, freqs);
    }


    public MicrosatelliteModel getSubModel(){
        return subModel;
    }

    public double[] getRates(){
        return super.getRates();
    }

    public void computeStationaryDistribution(){
       super.computeStationaryDistribution();
    }

    public abstract void setupInfinitesimalRates();


    public double getLogOneTransitionProbabilityEntry(double distance, int parentState, int childState){
        return Math.log(getOneTransitionProbabilityEntry(distance, parentState, childState));
    }

    public double getOneTransitionProbabilityEntry(double distance, int parentState, int childState){
       double probability = 0.0;
       double temp;

        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        if (!wellConditioned) {
            return 0.0;
        }

        double [] iexp = new double[stateCount];
        for(int i = 0; i < stateCount; i++){
            if(EvalImag[i] == 0){                
                temp = Math.exp(distance*(Eval[i]));
                iexp[i] = temp*Ievc[i][childState];
            }else{
                throw new RuntimeException("imaginary eigen values");
            }
        }
        for(int i = 0; i < stateCount; i++){
            probability += Evec[parentState][i]*iexp[i];
        }

        if(probability <= 0.0){
            probability = minProb;
        }
        return probability;
    }

   
}