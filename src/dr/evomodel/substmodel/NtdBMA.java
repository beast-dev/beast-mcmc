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
        relativeRates[1] = Math.exp(logKappa.getValue(0)+modelChoose.getValue(0)*logTN.getValue(0));
        //rate CT value
        relativeRates[4] = Math.exp(logKappa.getValue(0)-modelChoose.getValue(0)*logTN.getValue(0));
        //rate AC value
        relativeRates[0] = Math.exp(modelChoose.getValue(1)*logAC.getValue(0));
        //rate AT value
        relativeRates[2] = Math.exp(modelChoose.getValue(1)*logAT.getValue(0));
        //rate GC value
        relativeRates[3] = Math.exp(modelChoose.getValue(1)*logGC.getValue(0));
        //rate GT value
        relativeRates[5] = Math.exp(modelChoose.getValue(1)*logGT.getValue(0));

    }

}
