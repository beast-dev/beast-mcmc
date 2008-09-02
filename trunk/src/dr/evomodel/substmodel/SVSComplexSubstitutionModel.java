package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class SVSComplexSubstitutionModel extends ComplexSubstitutionModel {

	public SVSComplexSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel, Parameter parameter, Parameter indicators) {
		super(name, dataType, rootFreqModel, parameter);
		this.indicators = indicators;
		addParameter(indicators);
	}

	protected double[] getRates() {
		double[] rates = infinitesimalRates.getParameterValues();
		double[] indies = indicators.getParameterValues();
		final int dim = rates.length;
		for (int i = 0; i < dim; i++)
			rates[i] *= indies[i];
		return rates;
	}

	private Parameter indicators;

}
