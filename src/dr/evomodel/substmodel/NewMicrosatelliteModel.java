package dr.evomodel.substmodel;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.evolution.datatype.Microsatellite;

/**
 * @author Chieh-Hsi Wu
 * Implementation of models by Watkins (2007)
 */
public class NewMicrosatelliteModel extends AbstractSubstitutionModel {

    public NewMicrosatelliteModel(Microsatellite msat, FrequencyModel rootFreqModel){
        super("NewMicrosatelliteModel", msat, rootFreqModel);

    }

    protected void ratesChanged() {};
    protected void setupRelativeRates(){};
    protected void frequenciesChanged() {};

}
