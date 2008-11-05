/*
 * GeneralSubstitutionModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.newsubstmodel;

import dr.evolution.datatype.*;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * <b>A general model of sequence substitution</b>. A general reversible class for any
 * data type.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GeneralSubstitutionModel.java,v 1.37 2006/05/05 03:05:10 alexei Exp $
 */
public class GeneralSubstitutionModel extends BaseSubstitutionModel {

	public static final String GENERAL_SUBSTITUTION_MODEL = "generalSubstitutionModel";
	public static final String DATA_TYPE = "dataType";
	public static final String RATES = "rates";
	public static final String RELATIVE_TO = "relativeTo";
	public static final String FREQUENCIES = "frequencies";

	/**
	 * the rate which the others are set relative to
	 */
	protected int ratesRelativeTo;

	/**
	 * constructor
	 *
	 * @param dataType the data type
	 */
	public GeneralSubstitutionModel(DataType dataType, FrequencyModel freqModel, Parameter parameter, int relativeTo) {

		super(GENERAL_SUBSTITUTION_MODEL, dataType, freqModel,
                new DefaultEigenSystem(dataType.getStateCount()));

		ratesParameter = parameter;
		if (ratesParameter != null) {
			addParameter(ratesParameter);
			ratesParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, ratesParameter.getDimension()));
		}
		setRatesRelativeTo(relativeTo);
	}

	/**
	 * constructor
	 *
	 * @param dataType the data type
	 */
	protected GeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, int relativeTo) {

		super(name, dataType, freqModel,
                new DefaultEigenSystem(dataType.getStateCount()));

		setRatesRelativeTo(relativeTo);
	}

	protected void frequenciesChanged() {
		// Nothing to precalculate
	}

	protected void ratesChanged() {
		// Nothing to precalculate
	}

    protected void setupRelativeRates(double[] rates) {
		for (int i = 0; i < rates.length; i++) {
			if (i == ratesRelativeTo) {
				rates[i] = 1.0;
			} else if (i < ratesRelativeTo) {
				rates[i] = ratesParameter.getParameterValue(i);
			} else {
				rates[i] = ratesParameter.getParameterValue(i - 1);
			}
		}
	}

	/**
	 * set which rate the others are relative to
	 */
	public void setRatesRelativeTo(int ratesRelativeTo) {
		this.ratesRelativeTo = ratesRelativeTo;
	}

	// *****************************************************************
	// Interface Model
	// *****************************************************************


	protected void storeState() {
	} // nothing to do

	/**
	 * Restore the additional stored state
	 */
	protected void restoreState() {
		updateMatrix = true;
	}

	protected void acceptState() {
	} // nothing to do

	/**
	 * Parses an element from an DOM document into a DemographicModel. Recognises
	 * ConstantPopulation and ExponentialGrowth.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return GENERAL_SUBSTITUTION_MODEL;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Parameter ratesParameter;

			XMLObject cxo = (XMLObject) xo.getChild(FREQUENCIES);
			FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

			DataType dataType = freqModel.getDataType();

			cxo = (XMLObject) xo.getChild(RATES);

			int relativeTo = cxo.getIntegerAttribute(RELATIVE_TO) - 1;
			if (relativeTo < 0) throw new XMLParseException(RELATIVE_TO + " must be 1 or greater");

			ratesParameter = (Parameter) cxo.getChild(Parameter.class);

			int rateCount = ((dataType.getStateCount() - 1) * dataType.getStateCount()) / 2;

			if (ratesParameter == null) {

				if (rateCount == 1) {
					// simplest model for binary traits...
				} else {
					throw new XMLParseException("No rates parameter found in " + getParserName());
				}
			} else if (ratesParameter.getDimension() != rateCount - 1) {
				throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount - 1) + " dimensions.");
			}

			return new GeneralSubstitutionModel(dataType, freqModel, ratesParameter, relativeTo);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A general reversible model of sequence substitution for any data type.";
		}

		public Class getReturnType() {
			return SubstitutionModel.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(FREQUENCIES, FrequencyModel.class),
				new ElementRule(RATES,
						new XMLSyntaxRule[]{
								AttributeRule.newIntegerRule(RELATIVE_TO, false, "The index of the implicit rate (value 1.0) that all other rates are relative to. In DNA this is usually G<->T (6)"),
								new ElementRule(Parameter.class, true)}
				)
		};

	};

	protected Parameter ratesParameter = null;
}