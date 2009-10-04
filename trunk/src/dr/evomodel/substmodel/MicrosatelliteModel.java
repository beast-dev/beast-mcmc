package dr.evomodel.substmodel;


import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.inference.model.Model;


/**
 * @author Chieh-Hsi Wu
 *
 * An abstract class for microsatellite models
 */
public abstract class MicrosatelliteModel extends ComplexSubstitutionModel{
    protected OnePhaseModel subModel = null;
    protected boolean isNested = false;
    protected double[][] infinitesimalRateMatrix = null;



    /**
     * Constructor
     * @param name              Model name
     * @param msat              Microsatellite data type
     * @param rootFreqModel     Frequency model
     * @param parameter         Infinitesimal rates
     */
    public MicrosatelliteModel(String name, Microsatellite msat, FrequencyModel rootFreqModel, Parameter parameter) {
        super(name, msat, rootFreqModel, parameter);
        if(parameter == null){
            double[] q = new double[stateCount*(stateCount-1)];
            infinitesimalRates = new Parameter.Default(q);
        }

        infinitesimalRateMatrix = new double[stateCount][stateCount];
    }

    public void setToEqualFrequencies(){
        double[] freqs = new double[stateCount];
            for(int i = 0; i < freqs.length; i++){
                freqs[i] = 1.0/stateCount;
            }
        this.freqModel = new FrequencyModel(dataType, freqs);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        updateMatrix = true;
    }

    //store the infinitesimal rates in the vector to a matrix called amat
    public void storeIntoAmat(){
        amat = infinitesimalRateMatrix;
    }

    //matrix is already valid
    protected void makeValid(double[][] matrix, int dimension){}

    protected double getRate(int i, int j){
       return infinitesimalRateMatrix[i][j];
    }

    /*
     * Set up empirical frequencies
     */
    public void setupEmpiricalStationaryFrequencies(){
        int eigenValPos = 0;
        for(int i = 0; i < stateCount; i++){
            if(Eval[i] == 0){
                eigenValPos = i;
                break;
            }
        }
        double[] empFreq = new double[stateCount];
        for(int i = 0; i < stateCount; i++){
            empFreq[i] = Evec[i][eigenValPos]*Ievc[eigenValPos][i];

        }
        this.freqModel = new FrequencyModel(dataType, empFreq);
    }



    public double[] getRates(){
        return super.getRates();
    }

    public void computeStationaryDistribution(){
       super.computeStationaryDistribution();
    }

    public abstract void setupInfinitesimalRates();

    public Parameter getInfinitesimalRates(){
        return infinitesimalRates;
    }

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
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);
                iexp[i] = expatcosbt * Ievc[i][childState] + expatsinbt * Ievc[i2][childState];
                iexp[i2] = expatcosbt * Ievc[i2][childState] - expatsinbt * Ievc[i][childState];
                i ++;

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



    public MicrosatelliteModel getSubmodel(){
        return subModel;
    }

    public boolean isSubmodel(){
        return isNested;
    }

    public boolean hasSubmodel(){
        return subModel == null;
    }

}