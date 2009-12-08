package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * <b>A general model of sequence substitution with stochastic variable selection</b>. A general reversible class for any
 * data type.
 *
 * @author Marc Suchard
 * @version $Id: SVSIrreversibleSubstitutionModel.java,v 1.37 2006/05/05 03:05:10 msuchard Exp $
 */

public class SVSIrreversibleSubstitutionModel extends SVSGeneralSubstitutionModel {

	public SVSIrreversibleSubstitutionModel(DataType dataType, FrequencyModel freqModel,
	                                        FrequencyModel rootFreqModel, Parameter parameter,
	                                        Parameter indicator) {
		super(dataType, freqModel, parameter, indicator);

		this.rootFreq = rootFreqModel;
		addModel(rootFreq);
	}

	public FrequencyModel getFrequencyModel() {
		return rootFreq;
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// frequencyModel changed!

		// No need to update substitution model or tree likelihood when root frequencies change

		if (model == rootFreq)
			return;

		// If the stationary frequencies get changed, the treeLikelihood needs
		// to be recalculated, but treeLikelihood only listens to rootModel

		if (model == freqModel) {
			if (rootFreq != null)
				rootFreq.fireModelChanged();
			frequenciesChanged();
			updateMatrix = true;
		}
	}

	private FrequencyModel rootFreq;

    public boolean validState() {
        return true;
    }
    
}
