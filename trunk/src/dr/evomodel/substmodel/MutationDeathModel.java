/**
 * Package: MutationDeathModel
 * Description:
 *
 *
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 19, 2008
 * Time: 12:41:01 PM
 */

package dr.evomodel.substmodel;

import dr.evolution.datatype.MutationDeathType;
import dr.evolution.datatype.DataType;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class MutationDeathModel extends AbstractSubstitutionModel{

    public static final String MD_MODEL = "MutationDeathModel";
    private AbstractSubstitutionModel CTMCModel;
    private Parameter delParameter = null;
    protected double[] trMatrix;
    private Parameter baseSubModelFreq;
    private Parameter thisSubModelFreq;
    private Parameter mutationRate;



    MutationDeathModel(Parameter delParameter, DataType dT, AbstractSubstitutionModel evoModel,
                       FrequencyModel freqModel, Parameter mutationRate) {
        super(MD_MODEL, dT, freqModel);
        CTMCModel=evoModel;
        stateCount = freqModel.getFrequencyCount();
        this.delParameter=delParameter;
        this.dataType = dT;
        this.mutationRate = mutationRate;

        addParameter(delParameter);
        delParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        addParameter(mutationRate);
//        addModel(freqModel);
        if(evoModel!=null){
            addModel(evoModel.getFrequencyModel());
            addModel(evoModel);
        }
        trMatrix=new double[(stateCount-1)*(stateCount-1)];

        if(CTMCModel!=null){
            baseSubModelFreq = CTMCModel.getFrequencyModel().getFrequencyParameter();
        }
        else{
            baseSubModelFreq = new Parameter.Default(new double[]{1.0});
        }
        thisSubModelFreq = getFrequencyModel().getFrequencyParameter();

        double total = 0;
        for(int i=0; i<baseSubModelFreq.getDimension(); i++) {
            double value = thisSubModelFreq.getParameterValue(i);
            total += value;
            baseSubModelFreq.setParameterValue(i,value);
        }

        for(int i=0; i<baseSubModelFreq.getDimension(); i++) {
            baseSubModelFreq.setParameterValue(i,baseSubModelFreq.getParameterValue(i)/total);
        }

        thisSubModelFreq.setParameterValue(thisSubModelFreq.getDimension()-1,0.0);

        copyFrequencies();

        frequenciesChanged();
        ratesChanged();
        setupRelativeRates();
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

     protected void handleModelChangedEvent(Model model, Object object, int index) {

         if( object == baseSubModelFreq ) {
             copyFrequencies();
             fireModelChanged(object,index);
         } else if( model == CTMCModel ) {
             fireModelChanged(object,index);
         }
    }

    private void copyFrequencies() {
        for(int i=0; i<baseSubModelFreq.getDimension(); i++)
            thisSubModelFreq.setParameterValueQuietly(i,baseSubModelFreq.getParameterValue(i));

    }

    @Override
	protected void frequenciesChanged() {
	}

	@Override
	protected void ratesChanged() {
    }

	@Override
	protected void setupRelativeRates() {
    }
	
	public void getTransitionProbabilities(double distance, double[] matrix) {
        int i,j;
        // assuming that expected number of changes in CTMCModel is 1 per unit time
        // we are contributing s*deathRate number of changes per unit of time
        double deathProb=Math.exp(-distance*delParameter.getParameterValue(0));
        double mutationR=2*mutationRate.getParameterValue(0);
        double freqs[] = freqModel.getFrequencyParameter().getParameterValues();

        for(i=0;i<freqs.length-1;++i){
            mutationR*=freqs[i];
        }
        if(CTMCModel !=null){
            CTMCModel.getTransitionProbabilities(mutationR*distance,trMatrix);
        }
        else{
            trMatrix[0]=1.0;
        }

        for(i=0; i<stateCount-1;++i){
            for(j=0; j<stateCount-1;j++){
                matrix[i*(stateCount)+j]=trMatrix[i*(stateCount-1)+j]*deathProb;
            }
            matrix[i*(stateCount)+j]=(1.0-deathProb);
        }

        for(j=0; j<stateCount-1; ++j){
            matrix[stateCount*(stateCount-1)+j]=0.0;
        }

        matrix[stateCount*stateCount-1]=1.0;
/*        System.err.println("Transition matrix t = "+distance);
        for(i=0;i<stateCount;++i){
            System.err.print("[");
            for(j=0;j<stateCount; ++j){
                System.err.print(" "+matrix[i*stateCount+j]);
            }
            System.err.println(" ]");
        }*/
    }

    public static final String MUTATION_RATE ="mutationRate";
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

    public String getParserName() { return MD_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter dummyFreqParameter;
            Parameter delParam = (Parameter)xo.getChild(Parameter.class);

            Logger.getLogger("dr.evomodel").info("Creating MutationDeath substitution model.\n\tInitial death rate is "
                                                    + delParam.getParameterValue(0));

            MutationDeathType dT = (MutationDeathType) xo.getChild(MutationDeathType.class);

            AbstractSubstitutionModel evoModel = (AbstractSubstitutionModel)xo.getChild(AbstractSubstitutionModel.class);
            if(evoModel == null){  // Assuming pure survival model
                Logger.getLogger("dr.evomodel").info("\tSubstitutionModel not provided assuming pure death/survival model.");
                dummyFreqParameter = new Parameter.Default(new double[]{1.0,0.0});
            }
            else{
                dummyFreqParameter = new Parameter.Default(dT.getStateCount());
                double freqs[]=evoModel.getFrequencyModel().getFrequencies();
                for(int i = 0; i<freqs.length;++i){
                    dummyFreqParameter.setParameterValueQuietly(i,freqs[i]);
                }
                dummyFreqParameter.setParameterValueQuietly(dT.getStateCount()-1,0.0);
            }

            FrequencyModel dummyFrequencies = new FrequencyModel(dT,dummyFreqParameter);

            Parameter mutationRate;

            if(xo.hasSocket(MUTATION_RATE)){
                mutationRate = (Parameter) xo.getElementFirstChild(MUTATION_RATE);
            }
            else{
                mutationRate = new Parameter.Default(new double[]{1.0});
            }
            Logger.getLogger("dr.evomodel").info("\tInitial mutation rate is "+mutationRate.getParameterValue(0));

            return new MutationDeathModel(delParam, dT, evoModel,dummyFrequencies, mutationRate);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an instance of the MutationDeath model of CTMC evolution with deletions.";
        }

        public Class getReturnType() { return MutationDeathModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(AbstractSubstitutionModel.class,true),
            new ElementRule(Parameter.class),
            new ElementRule(MutationDeathType.class),
            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[] { new ElementRule(Parameter.class) },true)
        };

    };
}