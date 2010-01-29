package dr.evomodel.graphlikelihood;

import dr.evomodel.treelikelihood.GeneralLikelihoodCore;

/*
 * A class that supports calculating likelihoods on arbitrary ranges of sites
 * TODO: this should really be refactored into an abstract class that derives AbstractLikelihoodCore into AbstractGraphLikelihoodCore
 * and all the various likelihood cores should then derive from AbstractGraphLikelihoodCore.  Looks like total code volume could be
 * reduced by refactoring the matrixMap functions to merge them with the normal functions.
 */
public class GeneralGraphLikelihoodCore extends GeneralLikelihoodCore {

	public GeneralGraphLikelihoodCore(int stateCount) {
		super(stateCount);
	}

    /**
     * grows storage in partial likelihood arrays to accommodate more nodes.
     * reallocates and copies all storage
     * @param add           	the number of nodes to add
     */
    public void growNodeStorage(int add) {
        nodeCount += add;

        int[] currentMatricesIndicesTmp = new int[nodeCount];
        System.arraycopy(currentMatricesIndices, 0, currentMatricesIndicesTmp, 0, currentMatricesIndices.length);
        currentMatricesIndices = currentMatricesIndicesTmp;
        int[] storedMatricesIndicesTmp = new int[nodeCount];
        System.arraycopy(storedMatricesIndices, 0, storedMatricesIndicesTmp, 0, storedMatricesIndices.length);
        storedMatricesIndices = storedMatricesIndicesTmp;

        int[] currentPartialsIndicesTmp = new int[nodeCount];
        System.arraycopy(currentPartialsIndices, 0, currentPartialsIndicesTmp, 0, currentPartialsIndices.length);
        currentPartialsIndices = currentPartialsIndicesTmp;
        int[] storedPartialsIndicesTmp = new int[nodeCount];
        System.arraycopy(storedPartialsIndices, 0, storedPartialsIndicesTmp, 0, storedPartialsIndices.length);
        storedPartialsIndices = storedPartialsIndicesTmp;


        double [][][] partialsTmp = new double[partials.length][nodeCount][];
        int[][] statesTmp = new int[nodeCount][];
        double[][][] matricesTmp = new double[2][nodeCount][matrixCount * matrixSize];
    	double[][][] scalingFactorsTmp = new double[2][nodeCount][patternCount];
        for (int i = 0; i < nodeCount-add; i++) {
        	for(int j=0; j<2; j++)
        	{
        		partialsTmp[j][i] = partials[j][i];
                System.arraycopy(matrices[j][i], 0, matricesTmp[j][i], 0, matrixCount*matrixSize);
                if (useScaling) {
                    System.arraycopy(scalingFactors[j][i], 0, scalingFactorsTmp[j][i], 0, patternCount);
                }
        	}
            statesTmp[i] = states[i];
        }
        partials = partialsTmp;
        states = statesTmp;
        matrices = matricesTmp;

        if (useScaling) {
            scalingFactors = scalingFactorsTmp;
        }
        
        // assume the new nodes are internal nodes for now.
        for(int i=nodeCount-add; i<nodeCount; i++)
            createNodePartials(i);
        

    }

    //
    //
    // BEGIN OVERRIDES FROM AbstractLikelihoodCore
    //    
    
    /**
     * Calculates partial likelihoods at a node for a range of sites.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int left, int right) {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(
                        states[nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        states[nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], left, right);
            } else {
                calculateStatesPartialsPruning(states[nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex2]][nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], left, right);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(states[nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex1]][nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], left, right);
            } else {
                calculatePartialsPartialsPruning(partials[currentPartialsIndices[nodeIndex1]][nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex2]][nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], left, right);
            }
        }

        if (useScaling) {
        	// TODO: ensure that this works correctly when only part of the sites are updated!
            scalePartials(nodeIndex3);
        }
    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap  a map of which matrix to use for each pattern (can be null if integrating over categories)
     */
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap, int left, int right) {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(
                        states[nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        states[nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap, left, right);
            } else {
                calculateStatesPartialsPruning(
                        states[nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex2]][nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap, left, right);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(
                        states[nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex1]][nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap, left, right);
            } else {
                calculatePartialsPartialsPruning(
                        partials[currentPartialsIndices[nodeIndex1]][nodeIndex1], matrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndices[nodeIndex2]][nodeIndex2], matrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap, left, right);
            }
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    /**
     * Gets the partials for a particular node.
     *
     * @param nodeIndex   the node
     * @param outPartials an array into which the partials will go
     */
    public void getPartials(int nodeIndex, double[] outPartials, int left, int right) {
        double[] partials1 = partials[currentPartialsIndices[nodeIndex]][nodeIndex];

        System.arraycopy(partials1, left, outPartials, left, right-left);
    }

    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials, int left, int right) {
        calculateIntegratePartials(partials[currentPartialsIndices[nodeIndex]][nodeIndex], proportions, outPartials, left, right);
    }
    
    
    //
    //
    // BEGIN OVERRIDES FROM GeneralLikelihoodCore
    //    

    
    /**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
												int[] states2, double[] matrices2,
												double[] partials3, int left, int right)
	{
		int v = left;

		for (int l = 0; l < matrixCount; l++) {

			for (int k = left; k < right; k++) {

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
													double[] partials3, int left, int right)
	{

		double sum, tmp;

		int u = left;
		int v = left;

		for (int l = 0; l < matrixCount; l++) {
			for (int k = left; k < right; k++) {

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
													double[] partials3, int left, int right)
	{
		double sum1, sum2;

		int u = left;
		int v = left;

		for (int l = 0; l < matrixCount; l++) {

			for (int k = left; k < right; k++) {

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
												double[] partials3, int[] matrixMap, int left, int right)
	{
		int v = left;

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
													double[] partials3, int[] matrixMap, int left, int right)
	{

		double sum, tmp;

		int u = left;
		int v = left;

		for (int k = left; k < right; k++) {

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
													double[] partials3, int[] matrixMap, int left, int right)
	{
		double sum1, sum2;

		int u = left;
		int v = left;

		for (int k = left; k < right; k++) {

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
	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials, int left, int right)
	{

		int u = left;
		int v = left;
		for (int k = left; k < right; k++) {

			for (int i = 0; i < stateCount; i++) {

				outPartials[u] = inPartials[v] * proportions[0];
				u++;
				v++;
			}
		}


		for (int l = 1; l < matrixCount; l++) {
			u = left;

			for (int k = left; k < right; k++) {

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
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods, int left, int right)
	{
        int v = left;
		for (int k = left; k < right; k++) {

            double sum = 0.0;
			for (int i = 0; i < stateCount; i++) {

				sum += frequencies[i] * partials[v];
				v++;
			}
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}


}
