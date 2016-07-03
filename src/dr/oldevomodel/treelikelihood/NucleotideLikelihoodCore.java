/*
 * NucleotideLikelihoodCore.java
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

/**
 * NucleotideLikelihoodCore - An implementation of LikelihoodCore optimised
 * for nucleotides.
 *
 * @version $Id: NucleotideLikelihoodCore.java,v 1.12 2006/08/31 14:57:24 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */

@Deprecated // Switching to BEAGLE
public class NucleotideLikelihoodCore extends AbstractLikelihoodCore {

	/**
	 * Constructor
	 */
	public NucleotideLikelihoodCore() {

		super(4);

	}

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
												int[] states2, double[] matrices2,
												double[] partials3)
	{
		int v = 0;
		int u = 0;
		for (int j = 0; j < matrixCount; j++) {

			for (int k = 0; k < patternCount; k++) {

				int w = u;

				int state1 = states1[k];
				int state2 = states2[k];

				if (state1 < 4 && state2 < 4) {

					partials3[v] = matrices1[w + state1] * matrices2[w + state2];
					v++;	w += 4;
					partials3[v] = matrices1[w + state1] * matrices2[w + state2];
					v++;	w += 4;
					partials3[v] = matrices1[w + state1] * matrices2[w + state2];
					v++;	w += 4;
					partials3[v] = matrices1[w + state1] * matrices2[w + state2];
					v++;	w += 4;

				} else if (state1 < 4) {
					// child 2 has a gap or unknown state so don't use it

					partials3[v] = matrices1[w + state1];
					v++;	w += 4;
					partials3[v] = matrices1[w + state1];
					v++;	w += 4;
					partials3[v] = matrices1[w + state1];
					v++;	w += 4;
					partials3[v] = matrices1[w + state1];
					v++;	w += 4;

				} else if (state2 < 4) {
					// child 2 has a gap or unknown state so don't use it
					partials3[v] = matrices2[w + state2];
					v++;	w += 4;
					partials3[v] = matrices2[w + state2];
					v++;	w += 4;
					partials3[v] = matrices2[w + state2];
					v++;	w += 4;
					partials3[v] = matrices2[w + state2];
					v++;	w += 4;

				} else {
					// both children have a gap or unknown state so set partials to 1
					partials3[v] = 1.0;
					v++;
					partials3[v] = 1.0;
					v++;
					partials3[v] = 1.0;
					v++;
					partials3[v] = 1.0;
					v++;
				}
			}

			u += matrixSize;
		}
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3)
	{
		int u = 0;
		int v = 0;
		int w = 0;

		for (int l = 0; l < matrixCount; l++) {
			for (int k = 0; k < patternCount; k++) {

				int state1 = states1[k];

				if (state1 < 4) {

					double sum;

					sum =	matrices2[w] * partials2[v];
					sum +=	matrices2[w + 1] * partials2[v + 1];
					sum +=	matrices2[w + 2] * partials2[v + 2];
					sum +=	matrices2[w + 3] * partials2[v + 3];
					partials3[u] = matrices1[w + state1] * sum;	u++;

					sum =	matrices2[w + 4] * partials2[v];
					sum +=	matrices2[w + 5] * partials2[v + 1];
					sum +=	matrices2[w + 6] * partials2[v + 2];
					sum +=	matrices2[w + 7] * partials2[v + 3];
					partials3[u] = matrices1[w + 4 + state1] * sum;	u++;

					sum =	matrices2[w + 8] * partials2[v];
					sum +=	matrices2[w + 9] * partials2[v + 1];
					sum +=	matrices2[w + 10] * partials2[v + 2];
					sum +=	matrices2[w + 11] * partials2[v + 3];
					partials3[u] = matrices1[w + 8 + state1] * sum;	u++;

					sum =	matrices2[w + 12] * partials2[v];
					sum +=	matrices2[w + 13] * partials2[v + 1];
					sum +=	matrices2[w + 14] * partials2[v + 2];
					sum +=	matrices2[w + 15] * partials2[v + 3];
					partials3[u] = matrices1[w + 12 + state1] * sum;	u++;

					v += 4;

				} else {
					// Child 1 has a gap or unknown state so don't use it

					double sum;

					sum =	matrices2[w] * partials2[v];
					sum +=	matrices2[w + 1] * partials2[v + 1];
					sum +=	matrices2[w + 2] * partials2[v + 2];
					sum +=	matrices2[w + 3] * partials2[v + 3];
					partials3[u] = sum;	u++;

					sum =	matrices2[w + 4] * partials2[v];
					sum +=	matrices2[w + 5] * partials2[v + 1];
					sum +=	matrices2[w + 6] * partials2[v + 2];
					sum +=	matrices2[w + 7] * partials2[v + 3];
					partials3[u] = sum;	u++;

					sum =	matrices2[w + 8] * partials2[v];
					sum +=	matrices2[w + 9] * partials2[v + 1];
					sum +=	matrices2[w + 10] * partials2[v + 2];
					sum +=	matrices2[w + 11] * partials2[v + 3];
					partials3[u] = sum;	u++;

					sum =	matrices2[w + 12] * partials2[v];
					sum +=	matrices2[w + 13] * partials2[v + 1];
					sum +=	matrices2[w + 14] * partials2[v + 2];
					sum +=	matrices2[w + 15] * partials2[v + 3];
					partials3[u] = sum;	u++;

					v += 4;
				}
			}

			w += matrixSize;
		}
	}

	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3)
	{

		double sum1, sum2;

		int u = 0;
		int v = 0;
		int w = 0;

		for (int l = 0; l < matrixCount; l++) {
			for (int k = 0; k < patternCount; k++) {

				sum1 = matrices1[w] * partials1[v];
				sum2 = matrices2[w] * partials2[v];
				sum1 += matrices1[w + 1] * partials1[v + 1];
				sum2 += matrices2[w + 1] * partials2[v + 1];
				sum1 += matrices1[w + 2] * partials1[v + 2];
				sum2 += matrices2[w + 2] * partials2[v + 2];
				sum1 += matrices1[w + 3] * partials1[v + 3];
				sum2 += matrices2[w + 3] * partials2[v + 3];
				partials3[u] = sum1 * sum2; u++;

				sum1 = matrices1[w + 4] * partials1[v];
				sum2 = matrices2[w + 4] * partials2[v];
				sum1 += matrices1[w + 5] * partials1[v + 1];
				sum2 += matrices2[w + 5] * partials2[v + 1];
				sum1 += matrices1[w + 6] * partials1[v + 2];
				sum2 += matrices2[w + 6] * partials2[v + 2];
				sum1 += matrices1[w + 7] * partials1[v + 3];
				sum2 += matrices2[w + 7] * partials2[v + 3];
				partials3[u] = sum1 * sum2; u++;

				sum1 = matrices1[w + 8] * partials1[v];
				sum2 = matrices2[w + 8] * partials2[v];
				sum1 += matrices1[w + 9] * partials1[v + 1];
				sum2 += matrices2[w + 9] * partials2[v + 1];
				sum1 += matrices1[w + 10] * partials1[v + 2];
				sum2 += matrices2[w + 10] * partials2[v + 2];
				sum1 += matrices1[w + 11] * partials1[v + 3];
				sum2 += matrices2[w + 11] * partials2[v + 3];
				partials3[u] = sum1 * sum2; u++;

				sum1 = matrices1[w + 12] * partials1[v];
				sum2 = matrices2[w + 12] * partials2[v];
				sum1 += matrices1[w + 13] * partials1[v + 1];
				sum2 += matrices2[w + 13] * partials2[v + 1];
				sum1 += matrices1[w + 14] * partials1[v + 2];
				sum2 += matrices2[w + 14] * partials2[v + 2];
				sum1 += matrices1[w + 15] * partials1[v + 3];
				sum2 += matrices2[w + 15] * partials2[v + 3];
				partials3[u] = sum1 * sum2; u++;

				v += 4;
			}

			w += matrixSize;
		}

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

	/**
	 * Integrates partials across categories.
	 * @param inPartials the partials at the node to be integrated
	 * @param proportions the proportions of sites in each category
	 * @param outPartials an array into which the integrated partials will go
	 */
	public void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

		int u = 0;
		int v = 0;
		for (int k = 0; k < patternCount; k++) {

			outPartials[u] = inPartials[v] * proportions[0]; u++; v++;
			outPartials[u] = inPartials[v] * proportions[0]; u++; v++;
			outPartials[u] = inPartials[v] * proportions[0]; u++; v++;
			outPartials[u] = inPartials[v] * proportions[0]; u++; v++;
		}


		for (int j = 1; j < matrixCount; j++) {
			u = 0;
			for (int k = 0; k < patternCount; k++) {
				outPartials[u] += inPartials[v] * proportions[j]; u++; v++;
				outPartials[u] += inPartials[v] * proportions[j]; u++; v++;
				outPartials[u] += inPartials[v] * proportions[j]; u++; v++;
				outPartials[u] += inPartials[v] * proportions[j]; u++; v++;

			}
		}
	}

	/**
	 * Calculates site likelihoods at a node.
	 * @param partials the partials used to calculate the likelihoods
	 * @param frequencies an array of state frequencies
	 * @param outLogLikelihoods an array into which the likelihoods will go
	 */
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods)
	{
        int v = 0;
		for (int k = 0; k < patternCount; k++) {
			double sum = frequencies[0] * partials[v];	v++;
			sum += frequencies[1] * partials[v];	v++;
			sum += frequencies[2] * partials[v];	v++;
			sum += frequencies[3] * partials[v];	v++;
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}

}

