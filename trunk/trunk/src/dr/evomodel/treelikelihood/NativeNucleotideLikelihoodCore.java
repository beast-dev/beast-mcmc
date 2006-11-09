/*
 * NativeNucleotideLikelihoodCore.java
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

package dr.evomodel.treelikelihood;

/**
 * NativeNucleotideLikelihoodCore - An implementation of LikelihoodCore for nucleotides
 * that calls native methods for maximum speed. The native methods should be in
 * a shared library called "NucleotideLikelihoodCore" but the exact name will be system
 * dependent (i.e. "libNucleotideLikelihoodCore.so" or "NucleotideLikelihoodCore.dll").
 *
 * @version $Id: NativeNucleotideLikelihoodCore.java,v 1.14 2006/08/31 14:57:24 rambaut Exp $
 *
 * @author Andrew Rambaut
 */

public class NativeNucleotideLikelihoodCore extends AbstractLikelihoodCore {

	/**
	 * Constructor
	 */
	public NativeNucleotideLikelihoodCore() {

		super(4);

	}

	protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
												int[] states2, double[] matrices2,
												double[] partials3) {
        nativeStatesStatesPruning(states1, matrices1, states2, matrices2, patternCount, matrixCount, partials3);
    }

    protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3) {
        nativeStatesPartialsPruning(states1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3);
    }

    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3) {
        nativePartialsPartialsPruning(partials1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3);
    }

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
												int[] states2, double[] matrices2,
												double[] partials3, int[] matrixMap)
	{
		throw new RuntimeException("calculateStatesStatesPruning not implemented using matrixMap");
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3, int[] matrixMap)
	{
		throw new RuntimeException("calculateStatesStatesPruning not implemented using matrixMap");
	}

	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3, int[] matrixMap)
	{
		throw new RuntimeException("calculateStatesStatesPruning not implemented using matrixMap");
	}

	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials)
	{
		nativeIntegratePartials(inPartials, proportions, patternCount, matrixCount, outPartials);
	}


	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected native void nativeStatesStatesPruning(	int[] states1, double[] matrices1,
														int[] states2, double[] matrices2,
														int patternCount, int matrixCount,
														double[] partials3);
	protected native void nativeStatesPartialsPruning(	int[] states1, double[] matrices1,
														double[] partials2, double[] matrices2,
														int patternCount, int matrixCount,
														double[] partials3);
	protected native void nativePartialsPartialsPruning(double[] partials1, double[] matrices1,
														double[] partials2, double[] matrices2,
														int patternCount, int matrixCount,
														double[] partials3);
	public native void nativeIntegratePartials(			double[] partials, double[] proportions,
														int patternCount, int matrixCount,
														double[] outPartials);

	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods)
	{
        int v = 0;
		for (int k = 0; k < patternCount; k++) {
            double logScalingFactor = getLogScalingFactor(k);

			double sum = frequencies[0] * partials[v];	v++;
			sum += frequencies[1] * partials[v];	v++;
			sum += frequencies[2] * partials[v];	v++;
			sum += frequencies[3] * partials[v];	v++;
            outLogLikelihoods[k] = Math.log(sum) + logScalingFactor;
		}

        checkScaling();
    }

	public static boolean isAvailable() { return isNativeAvailable; }

	private static boolean isNativeAvailable = false;

	static {

		try {
			System.loadLibrary("NucleotideLikelihoodCore");

			// System.err.println("Native nucleotide likelihood core found");

			isNativeAvailable = true;
		} catch (UnsatisfiedLinkError e) {

			// System.err.println("Using Java nucleotide likelihood core");
		}

	}
}

