/*
 * DirichletProcessPriorLogger.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.app.bss.Utils;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

public class DirichletProcessPriorLogger implements Loggable {

//	private ParametricMultivariateDistributionModel baseModel;
	private Parameter precisionParameter;
	private Parameter categoriesParameter;
	private CompoundParameter uniquelyRealizedParameters;
	
	private int uniqueRealizationCount;
	private int realizationCount;

	private double[] categoryProbabilities;
	private int newCategoryIndex;
	private double meanForCategory;
	private double newX;

	public DirichletProcessPriorLogger(
			Parameter precisionParameter, //
			Parameter categoriesParameter, //
			CompoundParameter uniquelyRealizedParameters //
	) {

		this.precisionParameter = precisionParameter;
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;
		this.categoriesParameter = categoriesParameter;
		
		this.uniqueRealizationCount = uniquelyRealizedParameters.getDimension();
		this.realizationCount = categoriesParameter.getDimension();

	}// END: Constructor

	private double[] getCategoryProbs() {

		double[] probs = new double[uniqueRealizationCount];

		for (int i = 0; i < realizationCount; i++) {
			probs[(int) categoriesParameter.getParameterValue(i)]++;
		}// END: N loop

//		Utils.printArray(probs);
		
		for (int i = 0; i < uniqueRealizationCount; i++) {
			probs[i] = probs[i] / realizationCount;
		}// END: categoryCount loop

		// probs = new double[]{0.10, 0.10, 0.10, 0.10, 0.10,
		// 0.10, 0.10, 0.10, 0.10, 0.10};

		return probs;
	}// END: getCategoryProbs

	@Override
	public LogColumn[] getColumns() {

		// this.categoryProbabilities = getCategoryProbs();

		List<LogColumn> columns = new ArrayList<LogColumn>();

		columns.add(new NewLogger("x.new"));
		columns.add(new NewCategoryLogger("category.new"));
		columns.add(new NewMeanLogger("mean.new"));
		for (int i = 0; i < uniquelyRealizedParameters.getDimension(); i++) {
			columns.add(new ProbabilitiesLogger("pi.", i));
		}

		LogColumn[] rtnColumns = new LogColumn[columns.size()];

		return columns.toArray(rtnColumns);
	}// END: getColumns

	private void getNew() {

		this.categoryProbabilities = getCategoryProbs();

		this.newCategoryIndex = Utils.sample(categoryProbabilities);
		this.meanForCategory = uniquelyRealizedParameters
				.getParameterValue(newCategoryIndex);

		//TODO: generalize
		double sd = precisionParameter.getParameterValue(0);
		
//		System.out.println("FUBAR:" + sd);
		
		NormalDistribution nd = new NormalDistribution(meanForCategory, sd);
		this.newX = (Double) nd.nextRandom();

	}

	private class NewLogger extends NumberColumn {

		public NewLogger(String label) {
			super(label);
		}

		@Override
		public double getDoubleValue() {

			getNew();

			return newX;
		}// END: getDoubleValue

	}// END: NewLogger class

	private class NewMeanLogger extends NumberColumn {

		public NewMeanLogger(String label) {
			super(label);
		}

		@Override
		public double getDoubleValue() {

			return meanForCategory;
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

	private class ProbabilitiesLogger extends NumberColumn {

		private int i;

		public ProbabilitiesLogger(String label, int i) {
			super(label + i);
			this.i = i;
		}

		@Override
		public double getDoubleValue() {

			return categoryProbabilities[i];
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

	private class NewCategoryLogger extends NumberColumn {

		public NewCategoryLogger(String label) {
			super(label);
		}

		@Override
		public double getDoubleValue() {

			return newCategoryIndex;
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

}// END: class
