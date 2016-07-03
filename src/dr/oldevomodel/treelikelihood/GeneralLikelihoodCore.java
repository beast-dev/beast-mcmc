/*
 * GeneralLikelihoodCore.java
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
 * GeneralLikelihoodCore - An implementation of LikelihoodCore for any data
 *
 * @version $Id: GeneralLikelihoodCore.java,v 1.28 2006/08/31 14:57:24 rambaut Exp $
 *
 * @author Andrew Rambaut
 */

@Deprecated // Switching to BEAGLE
public class GeneralLikelihoodCore extends AbstractLikelihoodCore {

	/**
	 * Constructor
	 *
	 * @param stateCount number of states
	 */
	public GeneralLikelihoodCore(int stateCount) {
		super(stateCount);
	}

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
												int[] states2, double[] matrices2,
												double[] partials3)
	{
		int v = 0;

		for (int l = 0; l < matrixCount; l++) {

			for (int k = 0; k < patternCount; k++) {

				int state1 = states1[k];
				int state2 = states2[k];

				int w = l * matrixSize;

                if (state1 < stateCount && state2 < stateCount) {

					for (int i = 0; i < stateCount; i++) {

						partials3[v] = matrices1[w + state1] * matrices2[w + state2];

						v++;
						w += stateCount;
					}

				} else if (state1 < stateCount) {
					// child 2 has a gap or unknown state so treat it as unknown

					for (int i = 0; i < stateCount; i++) {

						partials3[v] = matrices1[w + state1];

						v++;
						w += stateCount;
					}
				} else if (state2 < stateCount) {
					// child 2 has a gap or unknown state so treat it as unknown

					for (int i = 0; i < stateCount; i++) {

						partials3[v] = matrices2[w + state2];

						v++;
						w += stateCount;
					}
				} else {
					// both children have a gap or unknown state so set partials to 1

					for (int j = 0; j < stateCount; j++) {
						partials3[v] = 1.0;
						v++;
					}
				}
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3)
	{

		double sum, tmp;

		int u = 0;
		int v = 0;

		for (int l = 0; l < matrixCount; l++) {
			for (int k = 0; k < patternCount; k++) {

				int state1 = states1[k];

                int w = l * matrixSize;

				if (state1 < stateCount) {


					for (int i = 0; i < stateCount; i++) {

						tmp = matrices1[w + state1];

						sum = 0.0;
						for (int j = 0; j < stateCount; j++) {
							sum += matrices2[w] * partials2[v + j];
							w++;
						}

						partials3[u] = tmp * sum;
						u++;
					}

					v += stateCount;
				} else {
					// Child 1 has a gap or unknown state so don't use it

					for (int i = 0; i < stateCount; i++) {

						sum = 0.0;
						for (int j = 0; j < stateCount; j++) {
							sum += matrices2[w] * partials2[v + j];
							w++;
						}

						partials3[u] = sum;
						u++;
					}

					v += stateCount;
				}
			}
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

		for (int l = 0; l < matrixCount; l++) {

			for (int k = 0; k < patternCount; k++) {

                int w = l * matrixSize;

				for (int i = 0; i < stateCount; i++) {

					sum1 = sum2 = 0.0;

					for (int j = 0; j < stateCount; j++) {
						sum1 += matrices1[w] * partials1[v + j];
						sum2 += matrices2[w] * partials2[v + j];

						w++;
					}

					partials3[u] = sum1 * sum2;
					u++;
				}
				v += stateCount;
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
												int[] states2, double[] matrices2,
												double[] partials3, int[] matrixMap)
	{
		int v = 0;

		for (int k = 0; k < patternCount; k++) {

			int state1 = states1[k];
			int state2 = states2[k];

			int w = matrixMap[k] * matrixSize;

			if (state1 < stateCount && state2 < stateCount) {

				for (int i = 0; i < stateCount; i++) {

					partials3[v] = matrices1[w + state1] * matrices2[w + state2];

					v++;
					w += stateCount;
				}

			} else if (state1 < stateCount) {
				// child 2 has a gap or unknown state so treat it as unknown

				for (int i = 0; i < stateCount; i++) {

					partials3[v] = matrices1[w + state1];

					v++;
					w += stateCount;
				}
			} else if (state2 < stateCount) {
				// child 2 has a gap or unknown state so treat it as unknown

				for (int i = 0; i < stateCount; i++) {

					partials3[v] = matrices2[w + state2];

					v++;
					w += stateCount;
				}
			} else {
				// both children have a gap or unknown state so set partials to 1

				for (int j = 0; j < stateCount; j++) {
					partials3[v] = 1.0;
					v++;
				}
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3, int[] matrixMap)
	{

		double sum, tmp;

		int u = 0;
		int v = 0;

		for (int k = 0; k < patternCount; k++) {

			int state1 = states1[k];

			int w = matrixMap[k] * matrixSize;

			if (state1 < stateCount) {

				for (int i = 0; i < stateCount; i++) {

					tmp = matrices1[w + state1];

					sum = 0.0;
					for (int j = 0; j < stateCount; j++) {
						sum += matrices2[w] * partials2[v + j];
						w++;
					}

					partials3[u] = tmp * sum;
					u++;
				}

				v += stateCount;
			} else {
				// Child 1 has a gap or unknown state so don't use it

				for (int i = 0; i < stateCount; i++) {

					sum = 0.0;
					for (int j = 0; j < stateCount; j++) {
						sum += matrices2[w] * partials2[v + j];
						w++;
					}

					partials3[u] = sum;
					u++;
				}

				v += stateCount;
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
													double[] partials2, double[] matrices2,
													double[] partials3, int[] matrixMap)
	{
		double sum1, sum2;

		int u = 0;
		int v = 0;

		for (int k = 0; k < patternCount; k++) {

			int w = matrixMap[k] * matrixSize;

			for (int i = 0; i < stateCount; i++) {

				sum1 = sum2 = 0.0;

				for (int j = 0; j < stateCount; j++) {
					sum1 += matrices1[w] * partials1[v + j];
					sum2 += matrices2[w] * partials2[v + j];

					w++;
				}

				partials3[u] = sum1 * sum2;
				u++;
			}
			v += stateCount;
		}
	}

	/**
	 * Integrates partials across categories.
     * @param inPartials the array of partials to be integrated
	 * @param proportions the proportions of sites in each category
	 * @param outPartials an array into which the partials will go
	 */
	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials)
	{

		int u = 0;
		int v = 0;
		for (int k = 0; k < patternCount; k++) {

			for (int i = 0; i < stateCount; i++) {

				outPartials[u] = inPartials[v] * proportions[0];
				u++;
				v++;
			}
		}


		for (int l = 1; l < matrixCount; l++) {
			u = 0;

			for (int k = 0; k < patternCount; k++) {

				for (int i = 0; i < stateCount; i++) {

					outPartials[u] += inPartials[v] * proportions[l];
					u++;
					v++;
				}
			}
		}
	}

	/**
	 * Calculates patten log likelihoods at a node.
	 * @param partials the partials used to calculate the likelihoods
	 * @param frequencies an array of state frequencies
	 * @param outLogLikelihoods an array into which the likelihoods will go
	 */
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods)
	{
        int v = 0;
		for (int k = 0; k < patternCount; k++) {

            double sum = 0.0;
			for (int i = 0; i < stateCount; i++) {

				sum += frequencies[i] * partials[v];
				v++;
			}
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}
}
