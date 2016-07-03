/*
 * NativeCovarionLikelihoodCore.java
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

package dr.oldevomodel.treelikelihood;

@Deprecated // Switching to BEAGLE
public class NativeCovarionLikelihoodCore extends AbstractLikelihoodCore {

	public NativeCovarionLikelihoodCore() {
		super(8);
	}

	protected void calculateIntegratePartials(double[] inPartials,
	                                          double[] proportions, double[] outPartials) {
		nativeIntegratePartials(inPartials, proportions, patternCount, matrixCount, outPartials);
	}

	protected void calculatePartialsPartialsPruning(double[] partials1,
	                                                double[] matrices1, double[] partials2, double[] matrices2,
	                                                double[] partials3) {
		nativePartialsPartialsPruning(partials1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3);
	}

	protected void calculateStatesPartialsPruning(int[] states1,
	                                              double[] matrices1, double[] partials2, double[] matrices2,
	                                              double[] partials3) {

		nativeStatesPartialsPruning(states1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3);
	}

	protected void calculateStatesStatesPruning(int[] states1,
	                                            double[] matrices1, int[] states2, double[] matrices2,
	                                            double[] partials3) {

		nativeStatesStatesPruning(states1, matrices1, states2, matrices2, patternCount, matrixCount, partials3);
	}

	protected void calculatePartialsPartialsPruning(double[] partials1,
	                                                double[] matrices1, double[] partials2, double[] matrices2,
	                                                double[] partials3, int[] matrixMap) {
		throw new RuntimeException("not implemented using matrixMap");
	}

	protected void calculateStatesStatesPruning(int[] states1,
	                                            double[] matrices1, int[] states2, double[] matrices2,
	                                            double[] partials3, int[] matrixMap) {
		throw new RuntimeException("not implemented using matrixMap");
	}

	protected void calculateStatesPartialsPruning(int[] states1,
	                                              double[] matrices1, double[] partials2, double[] matrices2,
	                                              double[] partials3, int[] matrixMap) {
		throw new RuntimeException("not implemented using matrixMap");
	}

	public void calculateLogLikelihoods(double[] partials,
	                                    double[] frequencies, double[] outLogLikelihoods) {

		int v = 0;
		for (int k = 0; k < patternCount; k++) {

			double sum = frequencies[0] * partials[v];
			v++;
			sum += frequencies[1] * partials[v];
			v++;
			sum += frequencies[2] * partials[v];
			v++;
			sum += frequencies[3] * partials[v];
			v++;

			sum += frequencies[4] * partials[v];
			v++;
			sum += frequencies[5] * partials[v];
			v++;
			sum += frequencies[6] * partials[v];
			v++;
			sum += frequencies[7] * partials[v];
			v++;

			outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}

	public native void nativeIntegratePartials(double[] partials, double[] proportions,
	                                           int patternCount, int matrixCount,
	                                           double[] outPartials);

	protected native void nativePartialsPartialsPruning(double[] partials1, double[] matrices1,
	                                                    double[] partials2, double[] matrices2,
	                                                    int patternCount, int matrixCount,
	                                                    double[] partials3);

	protected native void nativeStatesPartialsPruning(int[] states1, double[] matrices1,
	                                                  double[] partials2, double[] matrices2,
	                                                  int patternCount, int matrixCount,
	                                                  double[] partials3);

	protected native void nativeStatesStatesPruning(int[] states1, double[] matrices1,
	                                                int[] states2, double[] matrices2,
	                                                int patternCount, int matrixCount,
	                                                double[] partials3);

	public static boolean isAvailable() {
		return isNativeAvailable;
	}

	private static boolean isNativeAvailable = false;

	static {
		try {
			System.loadLibrary("CovarionLikelihoodCore");
			isNativeAvailable = true;
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Using Java general likelihood core " + e.toString());
		}
	}

}
