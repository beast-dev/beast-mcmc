package dr.evomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 * 
 * Implements the LinearBiasModel of microsatellites.
 */
public class LinearBiasModel extends OnePhaseModel{

    /*
     *Parameters for setting up the infinitesimal rate matrix.
     */
    private Parameter biasConst;
    private Parameter biasLin;
    private ArrayList<Parameter> submodelParameters = null;
    private boolean estimateSubmodelParams = false;
    private boolean updateSubmodelRates = false;
    private boolean inLogitSpace = false;

    public double delta = 1e-15;

    public static final String LINEAR_BIAS_MODEL = "LINEARBIASModel";


    /**
     * Constructor
     *
     * @param microsatellite            microsatellite data type
     * @param freqModel                 user specified initial equilibrium frequencies
     * @param submodel                  submodel of this linear bias model
     * @param biasConst                 bias constant parameter
     * @param biasLinear                bias linear parameter
     * @param inLogitSpace              indicates whether the bias parameters are in the logit space
     * @param estimateSubmodelParams    inidicate whether the parameters of submodel will be estimated
     * @param isSubmodel                inidicate whether this model is a submodel of another microsatellite model
     */
    public LinearBiasModel(
            Microsatellite microsatellite,
            FrequencyModel freqModel,
            OnePhaseModel submodel,
            Parameter biasConst,
            Parameter biasLinear,
            boolean inLogitSpace,
            boolean estimateSubmodelParams,
            boolean isSubmodel){

        super(LINEAR_BIAS_MODEL, microsatellite, freqModel, null);

        isNested = isSubmodel;
        this.subModel = submodel;
        this.estimateSubmodelParams = estimateSubmodelParams;
        if(this.estimateSubmodelParams){
            submodelParameters = new ArrayList<Parameter>();
            for(int i = 0; i < subModel.getNestedParameterCount(); i++){

                if(isNested){
                    addVariable(subModel.getNestedParameter(i));
                }

                addParam(subModel.getNestedParameter(i));
                submodelParameters.add(subModel.getNestedParameter(i));
            }
            updateSubmodelRates = true;
        }



        //The default setting of the parameters gives infinitesimal rates with no directional bias.
        if(biasConst != null){
            this.biasConst = biasConst;
        }else{
            this.biasConst = new Parameter.Default(0.5);
        }


        if(biasLinear != null){
            biasLin = biasLinear;
        }else{
            biasLin = new Parameter.Default(0.0);
            System.out.println("The bias linear parameter of LinearBiasModel is not provided.");
            System.out.println("Using default value: "+biasLin.getParameterValue(0));
        }
        addParam(this.biasConst);
        addParam(this.biasLin);

        this.inLogitSpace = inLogitSpace;
        System.out.println("In logit-space: "+this.inLogitSpace);

        setupInfinitesimalRates();

        if(freqModel == null){
            System.out.println("Creating LinearBiasModel: using empirical frequencies");
            computeStationaryDistribution();
        }else{
            this.freqModel = freqModel;

        }
        addModel(this.freqModel);

    }
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(submodelParameters !=null && submodelParameters.indexOf((Parameter)variable) != -1){
            updateSubmodelRates = true;
        }
        updateMatrix = true;

    }

    /**
     * Setting up the infinitesimal Rates
     * If rates of the submodel at allele i are given by
     * X -> X + 1 at rate g(i)
     * X -> X - 1 at rate h(i)
     * where g(i) and h(i) are functions of the parent microsatellite length,
     * the rates are defined by the following equations:
     * X -> X + 1 at rate g(i)*alpha
     * X -> X - 1 at rate h(i)*(1 - alpha)
     * Alpha is defined as alpha(biasConst, biasLin, i) = max{0, min{1, u – v(i - k)}}
     */
    public void setupInfinitesimalRates(){
        if(updateSubmodelRates){
            subModel.setupInfinitesimalRates();
            updateSubmodelRates = false;
        }

        double biasConst = this.biasConst.getParameterValue(0);
        double biasLin = this.biasLin.getParameterValue(0);
        double expansionProb = 0.5;
        double rowSum;
        for(int i = 0; i < stateCount;i++){
            rowSum = 0.0;
            expansionProb = computeExpansionProb(biasConst,biasLin,i);
            if(expansionProb < delta){
                System.out.println("changing expan prob from " + expansionProb+ " to " + delta);
                expansionProb = delta;

            }else if (expansionProb > (1.0-delta)){

                System.out.println("changing expan prob from " + expansionProb+ " to " + (1.0-delta));
                expansionProb = 1.0-delta;
            }
            if(i - 1 > -1){
                infinitesimalRateMatrix[i][i - 1] = subModel.getRate(i,i-1)*(1.0 - expansionProb);
                rowSum = rowSum+infinitesimalRateMatrix[i][i - 1];
            }
            if(i + 1 < stateCount){
                infinitesimalRateMatrix[i][i + 1] = subModel.getRate(i,i+1)*expansionProb;
                rowSum = rowSum + infinitesimalRateMatrix[i][i + 1];
            }

            infinitesimalRateMatrix[i][i] = rowSum*-1;

        }

    }

    public double computeExpansionProb(double biasConst, double biasLin, int length){
        double expanProb = 0.5;
        if(inLogitSpace){
            double numerator = Math.exp(biasConst+biasLin*length);
            expanProb = numerator/(1+numerator);
        }else{
            expanProb = biasConst+biasLin*length;
        }
        return  expanProb;

    }

    public void setupMatrix(){
        setupInfinitesimalRates();
        super.setupMatrix();
    }



    public Parameter getBiasConstant(){
        return biasConst;
    }

    public Parameter getBiasLinearPercent(){
        return biasLin;
    }

    public boolean isEstimatingSubmodelParams(){
        return estimateSubmodelParams;
    }

    public boolean isInLogitSpace(){
        return inLogitSpace;
    }

}
