/*
 * ARGUniformPrior.java
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

package dr.evomodel.arg.coalescent;


import dr.evomodel.arg.ARGModel;
import dr.xml.*;
import org.apache.commons.math.util.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

public class ARGUniformPrior extends ARGCoalescentLikelihood {

	public static final String ARG_UNIFORM_PRIOR = "argUniformPrior";
	public static final String INITIAL_CALCULATIONS = "initialCalculations";
	public static final int INITIAL_DEFAULT = 5;


	/**
	 * First index represents <code>(number of taxa - 3)</code>.  Second index
	 * represents <code>(number of reassortments)</code>.
	 */
	public static final double[][] logARGCoalescentCount = {
			//3 external taxa
			{1.0986122886681098, 4.276666119016055, 8.265650165580329,
					12.882968485504067, 18.022948777876554, 23.61061585799083,
					29.589035299467483, 35.9136593725818, 42.54885778516801,
					49.465568765313016, 56.63967062533004, 64.05083469593846,
					71.6817017141028, 79.51727811078219, 87.54448452157624},
			//4 external taxa
			{2.8903717578961645, 7.049254841255837, 11.72325121817378,
					16.888557143580204, 22.489284935076803, 28.475109395164285,
					34.804239049061785, 41.44233437710982, 48.360992220084476,
					55.53645069653148, 62.948589014628716, 70.58017373916618,
					78.41629051730933, 86.44391155847957, 94.65156231759461},
			//5 external taxa
			{5.192956850890211, 10.127430784020902, 15.386536070471918,
					21.031563309246433, 27.041414773420907, 33.38480151281157,
					40.03192406503731, 46.95658863790172, 54.13620464676437,
					61.55131458837549, 69.18508096220405, 77.02283620054189,
					85.0517120025277, 93.26034012092553, 101.63861177667185},
			//6 external taxa
			{7.90100705199242, 13.480736877978641, 19.249258325378765,
					25.32054548939157, 31.698755104859107, 38.36733349428991,
					45.30602732354885, 52.49523955146093, 59.917151710215464,
					67.55588200276209, 75.39734855785784, 83.42905687418532,
					91.63988547749693, 100.01989325102055, 108.56015364817871},
			//7 external taxa
			{10.945529489715843, 17.07892753271249, 23.29838718386135,
					29.753407176990198, 36.467041602426036, 43.43431706668038,
					50.64269194407989, 58.077999881992824, 65.7264083635128,
					73.57505740705963, 81.61221669682189, 89.8272615722345,
					98.21058447998666, 106.75348902699808, 115.44808566836372},
			//8 external taxa
			{14.277733999891046, 20.896472983408266, 27.5204956568849,
					34.32463786712059, 41.34646380715327, 48.590254195514014,
					56.04978852563933, 63.71543385947847, 71.57673212169905,
					79.62341044573522, 87.84577046007783, 96.23481579273648,
					104.78226604090887, 113.48052233027693, 122.32261419599219},
			//9 external taxa
			{17.861252938347157, 24.912242385415205, 31.903176098451674,
					39.02769936297709, 46.334764858303096, 53.836163974100344,
					61.53092519032464, 69.4132445398445, 77.47555198269822,
					85.70981250210274, 94.10810078400579, 102.66285174696894,
					111.36695890554311, 120.21379947792728, 129.1972243144967},
			//10 external taxa
			{21.667915428117478, 29.108649135506738, 36.43531279311537,
					43.855945609774885, 51.42860726365353, 59.171269494486445,
					67.08736948530374, 75.17436983012065, 83.4271948608094,
					91.83976176721744, 100.40570480277214, 109.11872661806457,
					117.97276669219758, 126.96207643311605, 136.08124559689557},
			//11 external taxa
			{25.67524861334995, 33.470895149684544, 41.10703544874183,
					48.802991754828916, 56.62422391351782, 64.59385162827832,
					72.71905241455046, 81.00010003398279, 89.43409011267886,
					98.0166551565224, 106.74281045596386, 115.60738869034444,
					124.6052669677749, 133.73148419256228, 142.98129861277982},
			//12 external taxa
			{29.864903355376374, 37.986383730127116, 45.90958634624165,
					53.86285360119818, 61.91774245806648, 70.10170706136994,
					78.42512661349667, 86.890712820353, 95.49746436335495,
					104.24253036074631, 113.12215615381358, 122.13218720871606,
					131.26834550980237, 140.52638287915215, 149.90216524590267},
			//13 external taxa
			{34.22161218206597, 42.644274889635966, 50.83517511899286,
					59.02998615127507, 67.30534956176908, 75.69240403974487,
					84.20429128794189, 92.84585184006139, 101.61776251108103,
					110.51852076251443, 119.54547309663305, 128.69537701990959,
					137.96471924891867, 147.34989978600277, 156.84733922685047},
			//14 external taxa
			{38.732471688582805, 47.43514810013058, 55.87684474986188,
					64.2992772236583, 72.78337269607047, 81.36342882144557,
					90.05498796388106, 98.86476107289788, 107.79491286633757,
					116.84514548111035, 126.01379567653075, 135.29844575074105,
					144.69627736969085, 154.20428208422834, 163.8193886161766},
			//15 external taxa
			{43.38643203874033, 52.35074398686485, 61.028354989069626,
					69.6660218145462, 78.34831802900398, 87.11227051188368,
					95.97552137406072, 104.94643415554035, 114.02849935687857,
					123.22250024459373, 132.52766774855348, 141.9423328912873,
					151.46431112672838, 161.0911356940919, 170.82020228928107}
	};

	private ArrayList<Double> argNumber;

	public ARGUniformPrior(ARGModel arg, int max, int initial) {
		super(ARG_UNIFORM_PRIOR, arg, max);

		addModel(arg);

		argNumber = new ArrayList<Double>(15);

		if (arg.getExternalNodeCount() - 3 < logARGCoalescentCount.length) {
			Logger.getLogger("dr.evomodel").info("Creating ARGUniformPrior using stored arg counts");
			for (int i = 0, n = arg.getExternalNodeCount() - 3; i < logARGCoalescentCount[n].length; i++)
				argNumber.add(logARGCoalescentCount[n][i]);
		} else {
			Logger.getLogger("dr.evomodel").info("Creating ARGUniformPrior by calculating arg counts");
			for (int i = 0, n = arg.getExternalNodeCount(); i < initial; i++) {
				argNumber.add(logNumberARGS(n, i));
			}
		}
	}

	public double getLogARGNumber(int i) {
		if (i >= argNumber.size()) {
			argNumber.add(logNumberARGS(arg.getExternalNodeCount(), i));
		}
		return argNumber.get(i);
	}

	public double getLogLikelihood() {
		if (likelihoodKnown) {
			return logLikelihood;
		}

		likelihoodKnown = true;
		logLikelihood = calculateLogLikelihood();


		if (arg.getReassortmentNodeCount() > maxReassortments)
			logLikelihood = Double.NEGATIVE_INFINITY;
		else
			logLikelihood = calculateLogLikelihood();

		if (!currentARGValid(true)) {
			logLikelihood = Double.NEGATIVE_INFINITY;
		}

		return logLikelihood;
	}

	public double calculateLogLikelihood() {


		double treeHeight = arg.getNodeHeight(arg.getRoot());
		int internalNodes = arg.getInternalNodeCount() - 1;


		double logLike = logFactorial(internalNodes) - (double) internalNodes * Math.log(treeHeight)
				- getLogARGNumber(arg.getReassortmentNodeCount());

		assert !Double.isInfinite(logLike) && !Double.isNaN(logLike);


		return logLike;
	}

	private double logFactorial(int n) {
		double rValue = 0;

		for (int i = n; i > 0; i--) {
			rValue += Math.log(i);
		}
		return rValue;
	}


	private int numberARGS(int taxa, int argNumber) {
		int x = taxa;
		int n = 2 * argNumber + taxa - 1;

		return shurikoRecursion(x, n);
	}

	private int shurikoRecursion(int x, int n) {
		int a = 0;
		if (x == 0) {
			a = 0;
		} else if (x == 1) {
			if (n == 0) {
				a = 1;
			} else {
				a = 0;
			}
		} else if (n == 0) {
			if (x == 1) {
				a = 1;
			} else {
				a = 0;
			}
		} else if (x == n + 1) {
			a = x * (x - 1) / 2 * shurikoRecursion(x - 1, n - 1);
		} else {
			a = x * shurikoRecursion(x + 1, n - 1) + x * (x - 1) / 2 * shurikoRecursion(x - 1, n - 1);
		}
		return a;
	}

	public static double logNumberARGS(int start, int reassortments) {
		Logger.getLogger("dr.evomodel").warning("Calculating ARG count for " + reassortments  + " reassortments.  This may take awhile");
		
		if (reassortments == 0) {
			double a = 0;
			for (int i = start; i > 2; i--) {
				a += Math.log(i * (i - 1) / 2.0);
			}
			return a;
		}

		int[] max = new int[start - 3 + reassortments * 2];
		int[] x = new int[max.length];

		int i = 0;
		while (i < reassortments) {
			x[i] = max[i] = 1;
			i++;
		}
		while (i < max.length) {
			x[i] = max[i] = -1;
			i++;
		}
		double before = 100;

		double approx = 0;
		while (x[0] != -9 && !stopCombination(x, start)) {
			if (testCombination(x, start)) {
				before = approx;
				int[] y = generateValues(x, start);
				approx += reduceThenDivide(y, generateValues(max, start));
				before = approx - before;
			}
			nextCombination(x);
		}

		approx = Math.log(approx);

		int[] y = new int[max.length + 2];
		for (i = 0; i < max.length; i++)
			y[i] = max[i];
		y[y.length - 2] = y[y.length - 1] = -1;

		max = generateValues(y, start);

		for (int k = 0; k < y.length; k++)
			approx += Math.log(max[k]);

		return approx;

	}

	private static double reduceThenDivide(int[] top, int[] bottom) {

		if (false) {
			for (int i = 0; i < top.length; i++) {
				for (int j = 0; j < bottom.length; j++) {
					int gcd = MathUtils.gcd(top[i], bottom[j]);

					if (gcd > 1) {
						top[i] = top[i] / gcd;
						bottom[j] = bottom[j] / gcd;
					}
				}

			}

		}

		Arrays.sort(top);
		Arrays.sort(bottom);


		double a = 1;
		for (int i = 0; i < top.length; i++)
			a *= (double) top[i] / bottom[i];
		return a;

	}

	private static int[] generateValues(int[] x, int start) {
		int[] y = new int[x.length];

		for (int i = 0; i < x.length; i++) {
			if (x[i] == 1)
				y[i] = start;
			else
				y[i] = start * (start - 1) / 2;

			start += x[i];

		}
		return y;
	}

	private static boolean testCombination(int[] x, int start) {

		for (int i = 0; i < x.length; i++) {
			start += x[i];
			if (start == 1)
				return false;
		}
		return true;

	}

	private static boolean stopCombination(int[] x, int start) {
		for (int i = 0; i < x.length; i++) {
			if (x[i] == -1) {
				start--;
				if (start == 1) {
					return true;
				}
			} else {
				break;
			}
		}
		return false;
	}

	private static void nextCombination(int[] x) {
		if (x[x.length - 1] == -1) {
			int i = x.length - 1;

			while (i > -1) {
				if (x[i] == 1) {
					x[i] = -1;
					x[i + 1] = 1;
					return;
				} else
					i--;

			}
		} else {
			int endOnes = 0;
			int i = x.length - 1;
			while (x[i] == 1) {
				endOnes++;
				i--;
			}
			int nextOne = -1;
			while (i > -1) {
				if (x[i] == 1) {
					nextOne = i;
					break;
				} else
					i--;

			}
			if (nextOne == -1) {
				x[0] = -9;
				return;
			}

			x[nextOne] = -1;
			x[nextOne + 1] = 1;

			for (i = 0; i < endOnes; i++)
				x[i + nextOne + 2] = 1;

			i = nextOne + 2 + endOnes;

			while (i < x.length) {
				x[i] = -1;
				i++;
			}
		}
	}


	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserDescription() {
			return "A uniform prior for an ARG model";
		}

		public Class getReturnType() {
			return ARGUniformPrior.class;
		}

		public String getParserName() {
			return ARG_UNIFORM_PRIOR;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(ARGModel.class),

		};

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			ARGModel argModel = (ARGModel) xo.getChild(ARGModel.class);

			int max = Integer.MAX_VALUE;
			if (xo.hasAttribute(MAX_REASSORTMENTS)) {
				max = xo.getIntegerAttribute(MAX_REASSORTMENTS);
			}

			int initial = INITIAL_DEFAULT;
			if (xo.hasAttribute(INITIAL_CALCULATIONS)) {
				initial = xo.getIntegerAttribute(INITIAL_CALCULATIONS);
			}

			return new ARGUniformPrior(argModel, max, initial);
		}


	};


}
