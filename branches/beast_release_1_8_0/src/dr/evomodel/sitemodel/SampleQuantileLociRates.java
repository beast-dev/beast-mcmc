package dr.evomodel.sitemodel;

import dr.inference.model.*;
import dr.inference.distribution.ParametricDistributionModel;
import dr.evomodelxml.sitemodel.SampleQuantileLociRatesParser;

/**
 * @author Chieh-Hsi Wu
 * Sampling relative loci rates from an underlying distribution by sampling quantiles.
 */
public class SampleQuantileLociRates extends AbstractModel {
    private CompoundParameter lociRates;
    private Parameter rateQuantileParameter;
    private ParametricDistributionModel distrModel;
    private double normalizeRateTo;
    private boolean normalize;
    private double scaleFactor;
    private double[] rates;



    public SampleQuantileLociRates(
            CompoundParameter lociRates,
            Parameter rateQuantileParameter,
            ParametricDistributionModel model,
            boolean normalize,
            double normalizeLociRateTo) {
        super(SampleQuantileLociRatesParser.SAMPLE_QUANTILE_LOCI_RATES);

        this.lociRates = lociRates;

        this.rateQuantileParameter = rateQuantileParameter;
        //Force the boundaries of rateCategoryParameter to match the category count
        //Parameter.DefaultBounds bound = new Parameter.DefaultBounds(0.99, 0.01, rateQuantileParameter.getDimension());
        //this.rateQuantileParameter.addBounds(bound);

        this.distrModel = model;
        this.normalizeRateTo = normalizeLociRateTo;
        this.normalize = normalize;
        rates = new double[rateQuantileParameter.getDimension()];
        setupRates();
        setupRelativeRates();

        addModel(distrModel);
        addVariable(this.rateQuantileParameter);
    }

    private void setupRates(){
        for(int i = 0; i < rates.length; i++){
            rates[i]= distrModel.quantile(rateQuantileParameter.getParameterValue(i));
            //System.err.println("setupRates: "+rates[i]);
        }
    }

    private void setupRates(int rateQuantileParameterIndex){
        rates[rateQuantileParameterIndex]= distrModel.quantile(rateQuantileParameter.getParameterValue(rateQuantileParameterIndex));
    }
    
    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distrModel) {
            setupRates();
            setupRelativeRates();
            //System.out.println("speed investigation 1");
            fireModelChanged();
        }else if (model == rateQuantileParameter) {
            //System.out.println("speed investigation 2");
            setupRates(index);
            fireModelChanged(null, index);
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        //System.out.println("speed investigation 3");
        setupRates(index);
        setupRelativeRates();
        fireModelChanged(null, index);
    }

    protected void storeState() {
    }
    protected void acceptState() {
    }
    protected void restoreState() {
        //setupRates();
    }
    private void setupRelativeRates(){
        //System.err.println("computeFactor: "+normalize);
        if(normalize){
           computeFactor();
        }
        int lociCount = rateQuantileParameter.getDimension();
        for(int i = 0; i < lociCount; i ++){
            //System.err.println(rates[i]*scaleFactor);
            lociRates.setParameterValue(i,rates[i]*scaleFactor);
        }

    }

    private void computeFactor(){
        double sumRates = 0.0;
        int lociCount = rateQuantileParameter.getDimension();
        for(int i = 0; i < lociCount; i++){
            //System.err.println("computeFactor: "+rates[i]);
            sumRates += rates[i];
        }
        scaleFactor = normalizeRateTo/(sumRates/lociCount);
    }


}
