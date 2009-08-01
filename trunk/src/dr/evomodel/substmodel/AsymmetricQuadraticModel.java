package dr.evomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.app.beauti.options.MicrosatelliteModelType;

/**
 * @author Chieh-Hsi Wu
 *
 * Implements the Asymmetric Quadratic Model
 */
public class AsymmetricQuadraticModel extends OnePhaseModel{

    /*
     *Parameters for setting up the infinitesimal rate matrix.
     */
    private Parameter expanConst;
    private Parameter expanLin;
    private Parameter expanQuad;
    private Parameter contractConst;
    private Parameter contractLin;
    private Parameter contractQuad;



    /**
     * Constructor
     *
     * @param microsatellite    Mirosatellite data type
     * @param freqModel         Frequency model
     */
    public AsymmetricQuadraticModel(Microsatellite microsatellite, FrequencyModel freqModel){
        this(microsatellite, freqModel,
            null, null, null, null, null, null);

    }

    /**
     * Constructor
     *
     * @param microsatellite    Mirosatellite data type
     * @param freqModel         Frequency model
     * @param expanConst        Expansion constant
     * @param expanLin          Expansion linear coefficient
     * @param expanQuad         Expansion quadratic coefficient
     * @param contractConst     Contraction constant
     * @param contractLin       Contraction linear coefficient
     * @param contractQuad      Contraction quadratic coefficient
     */
    public AsymmetricQuadraticModel(Microsatellite microsatellite, FrequencyModel freqModel,
                    Parameter expanConst, Parameter expanLin, Parameter expanQuad,
                    Parameter contractConst, Parameter contractLin, Parameter contractQuad){
        
        super(MicrosatelliteModelType.ASYMQUAD.getXMLName(), microsatellite, freqModel,null);


        //The default setting of the parameters gives the same infinitesimal rates
        // as the truncated stepwise mutational model.
        this.expanConst = overrideDefault(new Parameter.Default(1.0), expanConst);
        this.expanLin = overrideDefault(new Parameter.Default(0.0), expanLin);
        this.expanQuad = overrideDefault(new Parameter.Default(0.0), expanQuad);
        this.contractConst = overrideDefault(new Parameter.Default(1.0), contractConst);
        this.contractLin = overrideDefault(new Parameter.Default(0.0), contractLin);
        this.contractQuad = overrideDefault(new Parameter.Default(0.0), contractQuad);

        addVariable(this.expanConst);
        addVariable(this.expanLin);
        addVariable(this.expanQuad);
        addVariable(this.contractConst);
        addVariable(this.contractLin);
        addVariable(this.contractQuad);

        setupInfinitesimalRates();

        //calculate the default frequencies when not provieded by the user.
        if(freqModel == null){
            System.out.println("Creating AysmmetricQuadraticModel: using empirical frequencies");
            computeStationaryDistribution();
        }else{
            this.freqModel = freqModel;
            addModel(this.freqModel);
        }
    }

    /*
     *  This method will override the default value of the parameter using the value specified by the user.
     */
    private Parameter overrideDefault(Parameter defaultParam, Parameter providedParam){
        if(providedParam != null)
            defaultParam = providedParam;
        return defaultParam;
    }


    /*
     * Setting up the infinitesimal Rates
     * The rates are defined by the following equations:
     * X -> X + 1 at rate u0 + u1(X - k) + u2(X - k)^2
     * X -> X - 1 at rate d0 + d1(X - k) + d2(X - k)^2
     */
    public void setupInfinitesimalRates(){

        double u0 = expanConst.getParameterValue(0);
        double u1 = expanLin.getParameterValue(0);
        double u2 = expanQuad.getParameterValue(0);

        double d0 = contractConst.getParameterValue(0);
        double d1 = contractLin.getParameterValue(0);
        double d2 = contractQuad.getParameterValue(0);

        double[] q = new double[stateCount*(stateCount-1)];
        int k = 0;
        for(int i = 0; i < stateCount; i++){
            for(int j = 0; j < stateCount; j++){
                if(i!=j){
                    if(j == i+1){

                        q[k] = u0+u1*i+u2*i*i;

                    }else if (j == i-1){

                        q[k] = d0+d1*i+d2*i*i;

                    }else{

                        q[k] = 0.0;

                    }
                    k++;
                }

            }
        }

        infinitesimalRates = new Parameter.Default(q);

    }
    

    public void setupMatrix(){
        setupInfinitesimalRates();
        computeStationaryDistribution();
        super.setupMatrix();
    }

    public Parameter getExpansionConstant(){
        return expanConst;
    }

    public Parameter getExpansionLinear(){
        return expanLin;
    }

    public Parameter getExpansionQuad(){
        return expanQuad;
    }

    public Parameter getContractionConstant(){
        return contractConst;
    }

    public Parameter getContractionLinear(){
        return contractLin;
    }

    public Parameter getContractionQuad(){
        return contractQuad;
    }
        
}
