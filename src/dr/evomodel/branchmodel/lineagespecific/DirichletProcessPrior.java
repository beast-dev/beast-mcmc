/*
 * DirichletProcessPrior.java
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
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

@SuppressWarnings("serial")
public class DirichletProcessPrior  extends AbstractModelLikelihood  {

	private static boolean VERBOSE = false;
	
	private Parameter categoriesParameter;
	private CompoundParameter uniquelyRealizedParameters;
	public ParametricMultivariateDistributionModel baseModel;
	private Parameter gamma;

	private int categoryCount;
	private int N;
	
	private boolean likelihoodKnown = false;
	private double logLikelihood;
	private final List<Double> cachedLogFactorials;

	public DirichletProcessPrior(Parameter categoriesParameter, //
			CompoundParameter uniquelyRealizedParameters, //
			ParametricMultivariateDistributionModel baseModel, //
			Parameter gamma //
	) {

		super("");

		// vector z of cluster assignments
		this.categoriesParameter = categoriesParameter;
		this.baseModel = baseModel;
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;
		this.gamma = gamma;
		
		// K clusters
		this.categoryCount = uniquelyRealizedParameters.getDimension();
//		this.categoryCount=Utils.findMaximum(categoriesParameter.getParameterValues()) + 1;
        this.N = categoriesParameter.getDimension();
		
		cachedLogFactorials = new ArrayList<Double>();
		cachedLogFactorials.add(0, 0.0);

		// add all
		this.addVariable(this.categoriesParameter);
		this.addVariable(this.gamma);
		
		this.addVariable(this.uniquelyRealizedParameters);
		
		if(baseModel != null) {
		this.addModel(baseModel);
		}
		
		this.likelihoodKnown = false;
	}// END: Constructor

	private double getLogFactorial(int i) {

		if ( cachedLogFactorials.size() <= i) {

			for (int j = cachedLogFactorials.size() - 1; j <= i; j++) {

				double logfactorial = cachedLogFactorials.get(j) + Math.log(j + 1);
				cachedLogFactorials.add(logfactorial);
			}

		}

		return cachedLogFactorials.get(i);
	}

	/**
	 * Assumes mappings start from index 0
	 * */
	private int[] getCounts() {

		// eta_k parameters (number of assignments to each category)
		int[] counts = new int[categoryCount];
		for (int i = 0; i < N; i++) {

			int category = getMapping(i);
			counts[category]++;

		}// END: i loop

		return counts;
	}// END: getCounts

	public double getGamma() {
		return gamma.getParameterValue(0);
	}
	
	private int getMapping(int i) {
		return (int) categoriesParameter.getParameterValue(i);
	}

	public double getLogDensity(Parameter parameter) {
		double value[] = parameter.getAttributeValue();
		return baseModel.logPdf(value);
	}

	public double getRealizedValuesLogDensity() {

		double total = 0.0;

		for (int i = 0; i < categoryCount; i++) {

			Parameter param = uniquelyRealizedParameters.getParameter(i);
			total += getLogDensity(param);
			
		}
		
		return  total;
	}//END: getRealizedValuesLogDensity

	public double getCategoriesLogDensity() {

		int[] counts = getCounts();
		
		if (VERBOSE) {
			Utils.printArray(counts);
		}
		
		double loglike = categoryCount * Math.log(getGamma());

		for (int k = 0; k < categoryCount; k++) {

			int eta = counts[k];

			if (eta > 0) {
				loglike += getLogFactorial(eta - 1);
			}

		}// END: k loop

		for (int i = 1; i <= N; i++) {
			loglike -= Math.log(getGamma() + i - 1);
		}// END: i loop

		return loglike;
	}// END: getPriorLoglike

	@Override
	public Model getModel() {
		return this;
	}

	@Override
	public double getLogLikelihood() {

		this.fireModelChanged();
		likelihoodKnown = false;
		if (!likelihoodKnown) {
			
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}

		return logLikelihood;
	}

	private double calculateLogLikelihood() {
//		getCounts();
		double loglike = getCategoriesLogDensity() + getRealizedValuesLogDensity();
	
		//TODOs
//		System.out.println(loglike);
		
		return loglike;
	}//END: calculateLogLikelihood

	@Override
	public void makeDirty() {
//		likelihoodKnown = false;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		likelihoodKnown = false;
	}

	public int getCategoryCount() {
		return categoryCount;
	}
	
	public Parameter getUniqueParameters() {
		return uniquelyRealizedParameters;
	}
	
	public Parameter getUniqueParameter(int index) {
		return uniquelyRealizedParameters.getParameter(index);
	}
	
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {

		if (variable == categoriesParameter) {
			
			this.fireModelChanged();
			
		} else if (variable == gamma) {

			// do nothing
			
			this.fireModelChanged();
			
		} else if (variable == uniquelyRealizedParameters) {

			likelihoodKnown = false;
			this.fireModelChanged();
			
		} else {

			throw new IllegalArgumentException("Unknown parameter");

		}

	}// END: handleVariableChangedEvent

	public void setVerbose() {
		VERBOSE = true;
	}
	
	@Override
	protected void storeState() {

	}

	@Override
	protected void restoreState() {
		likelihoodKnown = false;
	}

	@Override
	protected void acceptState() {

	}

	public static void main(String args[]) {

		testDirichletProcess(new double[] { 0, 1, 2 }, 3, 1.0, -Math.log(6.0));
		
		testDirichletProcess(new double[] { 0, 0, 1 }, 3, 1.0, -Math.log(6.0));
		
		testDirichletProcess(new double[] { 0, 1, 2, 3, 4 }, 5, 0.5, -6.851184927493743);
		
	}// END: main

	private static void testDirichletProcess(double[] mapping, int categoryCount,double gamma,
			double expectedLogL) {

		Parameter categoriesParameter = new Parameter.Default(mapping);
		Parameter gammaParameter = new Parameter.Default(gamma);

		CompoundParameter dummy = new CompoundParameter("dummy");

		for (int i = 0; i < categoryCount; i++) {
			dummy.addParameter(new Parameter.Default(1.0));
		}

		DirichletProcessPrior dpp = new DirichletProcessPrior(
				categoriesParameter, dummy, null, gammaParameter);
		dpp.setVerbose();
		
//		int[] counts = dpp.getCounts();
		System.out.println("lnL:          " + dpp.getCategoriesLogDensity());
		System.out.println("expected lnL: " + expectedLogL);
	}// END: testDirichletProcess
	
}// END: class

