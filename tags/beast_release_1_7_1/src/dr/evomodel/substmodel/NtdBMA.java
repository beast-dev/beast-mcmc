package dr.evomodel.substmodel;

import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.inference.model.Bounds;

/**
 * @author Chieh-Hsi Wu
 *
 * Model averaging of nucleotide substitution model
 *
 */
public class NtdBMA extends AbstractNucleotideModel{

    private Variable<Double> logKappa = null;
    private Variable<Double> logTN = null;
    private Variable<Double> logAC = null;
    private Variable<Double> logAT = null;
    private Variable<Double> logGC = null;
    private Variable<Double> logGT = null;
    private Variable<Integer> modelChoose = null;

    public static final int TN_INDEX = 0;
    public static final int GTR_INDEX = 1;

    public static final int ABSENT = 0;
    public static final int PRESENT = 1;

    public NtdBMA(
            Variable<Double> logKappa,
            Variable<Double> logTN,
            Variable<Double> logAC,
            Variable<Double> logAT,
            Variable<Double> logGC,
            Variable<Double> logGT,
            Variable<Integer> modelChoose,
            FrequencyModel freqModel){
        super("NucleotideBMA", freqModel);


        addVariable(logKappa);
        logKappa.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.logKappa = logKappa;

        addVariable(logTN);
        logTN.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.logTN = logTN;

        addVariable(logAC);
        logAC.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.logAC = logAC;

        addVariable(logAT);
        logAT.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.logAT = logAT;

        addVariable(logGC);
        logGC.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.logGC = logGC;

        addVariable(logGT);
        logGT.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        this.logGT = logGT;

        addVariable(modelChoose);
        modelChoose.addBounds(new Bounds.Int(modelChoose, 0, 1));

        
        this.modelChoose = modelChoose;
    }

    protected void setupRelativeRates() {

        //rate AG value
        relativeRates[1] = Math.exp(logKappa.getValue(0)+modelChoose.getValue(TN_INDEX)*logTN.getValue(0));
        //rate CT value
        relativeRates[4] = Math.exp(logKappa.getValue(0)-modelChoose.getValue(TN_INDEX)*logTN.getValue(0));
        //rate AC value
        relativeRates[0] = Math.exp(modelChoose.getValue(GTR_INDEX)*logAC.getValue(0));
        //rate AT value
        relativeRates[2] = Math.exp(modelChoose.getValue(GTR_INDEX)*logAT.getValue(0));
        //rate GC value
        relativeRates[3] = Math.exp(modelChoose.getValue(GTR_INDEX)*logGC.getValue(0));
        //rate GT value
        relativeRates[5] = Math.exp(modelChoose.getValue(GTR_INDEX)*logGT.getValue(0));

    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        boolean changeRateMatrix = false;
        
        if(variable == modelChoose){
            //changing the substitution model
            changeRateMatrix = true;
        }else if(variable == logKappa){
            //changing kappa value (which is in all models)
            changeRateMatrix = true;
        }else if(variable == logTN && modelChoose.getValue(TN_INDEX) == PRESENT){
            //changing the tn value when the current model is TN.
            changeRateMatrix = true;
        }else if(variable == logAC ||
                variable == logAT ||
                variable == logGC ||
                variable == logGT &&
                modelChoose.getValue(GTR_INDEX) == PRESENT && modelChoose.getValue(TN_INDEX) == PRESENT){
            //Changing any one of A<->C, A<->T, G<->C, G<->T rates while the current model is GTR.
            changeRateMatrix = true;
        }
        // relativeRates changed
        if(changeRateMatrix){
            super.handleVariableChangedEvent(variable, index, type);
        }
        changeRateMatrix = false;
    }

}
