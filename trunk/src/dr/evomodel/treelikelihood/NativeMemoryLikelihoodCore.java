package dr.evomodel.treelikelihood;

/*
 * NativeMemoryLikelihoodCore.java
 *
 * @author Marc Suchard
 * @author Dat Huynh
 *
 */

public class NativeMemoryLikelihoodCore implements LikelihoodCore {

	private static final long DEBUG = 100;

	protected int stateCount;
	protected int nodeCount;
	protected int patternCount;
	protected int partialsSize;
	protected int matrixSize;
	protected int matrixCount;

	protected boolean integrateCategories;

	protected long[][] ptrPartials;
	protected long[] ptrStates;
	protected long[][] ptrMatrices;

	protected int[] currentMatricesIndices;
	protected int[] storedMatricesIndices;
	protected int[] currentPartialsIndices;
	protected int[] storedPartialsIndices;

	protected boolean useScaling = false;

	protected double[][][] scalingFactors;

//    private double scalingThreshold = 1.0E-100;

	/**
	 * Constructor
	 *
	 * @param stateCount number of states
	 */
	public NativeMemoryLikelihoodCore(int stateCount) {
		this.stateCount = stateCount;

	}

	/**
	 * initializes partial likelihood arrays.
	 *
	 * @param nodeCount           the number of nodes in the tree
	 * @param patternCount        the number of patterns
	 * @param matrixCount         the number of matrices (i.e., number of categories)
	 * @param integrateCategories whether sites are being integrated over all matrices
	 */
	public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories) {

		this.nodeCount = nodeCount;
		this.patternCount = patternCount;
		this.matrixCount = matrixCount;

		this.integrateCategories = integrateCategories;

		if (integrateCategories) {
			partialsSize = patternCount * stateCount * matrixCount;
		} else {
			partialsSize = patternCount * stateCount;
		}

		ptrPartials = new long[2][nodeCount];

		currentMatricesIndices = new int[nodeCount];
		storedMatricesIndices = new int[nodeCount];

		currentPartialsIndices = new int[nodeCount];
		storedPartialsIndices = new int[nodeCount];

		ptrStates = new long[nodeCount];

		matrixSize = stateCount * stateCount;

		ptrMatrices = new long[2][nodeCount];

	}

	/**
	 * cleans up and deallocates arrays.
	 */
	public void finalize() throws Throwable {
		super.finalize();

		nodeCount = 0;
		patternCount = 0;
		matrixCount = 0;

		currentPartialsIndices = null;
		storedPartialsIndices = null;
		currentMatricesIndices = null;
		storedMatricesIndices = null;

//        scalingFactors = null;

		for (int i = 0; i < ptrPartials.length; i++) {
			for (int j = 0; j < ptrPartials[i].length; j++)
				wrapperFreeNativeMemoryArray(ptrPartials[i][j]);
		}

		for (int i = 0; i < ptrStates.length; i++)
			wrapperFreeNativeMemoryArray(ptrStates[i]);

		for (int i = 0; i < ptrMatrices.length; i++) {
			for (int j = 0; j < ptrMatrices[i].length; j++)
				wrapperFreeNativeMemoryArray(ptrMatrices[i][j]);
		}

	}


	public void printAllAddress() {

		System.err.println("Partials:");
		for (int i = 0; i < ptrPartials.length; i++) {
			for (int j = 0; j < ptrPartials[i].length; j++)
				System.err.println(ptrPartials[i][j] + " " + i + ":" + j);
		}

		System.err.println("States:");
		for (int i = 0; i < ptrStates.length; i++)
			System.err.println(ptrStates[i] + " " + i);

		System.err.println("Matrices:");
		for (int i = 0; i < ptrMatrices.length; i++) {
			for (int j = 0; j < ptrMatrices[i].length; j++)
				System.err.println(ptrMatrices[i][j] + " " + i + ":" + j);
		}


	}

	public void setUseScaling(boolean useScaling) {
		this.useScaling = useScaling;

		if (useScaling) {
			scalingFactors = new double[2][nodeCount][patternCount];
		}
	}

	/**
	 * Allocates partials for a node
	 */
	public void createNodePartials(int nodeIndex) {

		this.ptrPartials[0][nodeIndex] = wrapperAllocateNativeMemoryArray(partialsSize);
		this.ptrPartials[1][nodeIndex] = wrapperAllocateNativeMemoryArray(partialsSize);

	}


	/**
	 * Sets partials for a node
	 */
	public void setNodePartials(int nodeIndex, double[] partials) {

		if (this.ptrPartials[0][nodeIndex] == 0) {
			createNodePartials(nodeIndex);
		}

		double[] expandedPartials = partials;
		if (partials.length < partialsSize) {
			expandedPartials = new double[partials.length * matrixCount];
			int k = 0;
			for (int i = 0; i < matrixCount; i++) {
				System.arraycopy(partials, 0, expandedPartials, k, partials.length);
				k += partials.length;
			}
		}

		assert this.ptrPartials[0][nodeIndex] > DEBUG;

		setNativeMemoryArray(expandedPartials, 0, this.ptrPartials[0][nodeIndex], 0, expandedPartials.length);
	}


	/**
	 * Allocates states for a node
	 */
	public void createNodeStates(int nodeIndex) {

		this.ptrStates[nodeIndex] = wrapperAllocateNativeMemoryArray(patternCount);
	}

	/**
	 * Sets states for a node
	 */
	public void setNodeStates(int nodeIndex, int[] states) {

		if (this.ptrStates[nodeIndex] == 0) {
			createNodeStates(nodeIndex);
		}

		assert this.ptrStates[nodeIndex] > DEBUG;

		setNativeMemoryArray(states, 0, this.ptrStates[nodeIndex], 0, patternCount);
	}

	/**
	 * Gets states for a node
	 */
	public void getNodeStates(int nodeIndex, int[] states) {

		assert this.ptrStates[nodeIndex] > DEBUG;

		getNativeMemoryArray(this.ptrStates[nodeIndex], 0, states, 0, patternCount);
	}

	public void setNodeMatrixForUpdate(int nodeIndex) {
		currentMatricesIndices[nodeIndex] = 1 - currentMatricesIndices[nodeIndex];

	}


	/**
	 * Sets probability matrix for a node
	 */
	public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {

		if (ptrMatrices[currentMatricesIndices[nodeIndex]][nodeIndex] == 0) {
			ptrMatrices[currentMatricesIndices[nodeIndex]][nodeIndex] = wrapperAllocateNativeMemoryArray(matrixSize * matrixCount);
		}

		assert ptrMatrices[currentMatricesIndices[nodeIndex]][nodeIndex] > DEBUG;
		assert matrix != null;

		setNativeMemoryArray(matrix, 0, ptrMatrices[currentMatricesIndices[nodeIndex]][nodeIndex],
				matrixIndex * matrixSize, matrixSize);
	}

	/**
	 * Gets probability matrix for a node
	 */
	public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {

		assert ptrMatrices[currentMatricesIndices[nodeIndex]][nodeIndex] > DEBUG;

		assert matrix != null;

		getNativeMemoryArray(ptrMatrices[currentMatricesIndices[nodeIndex]][nodeIndex],
				matrixIndex * matrixSize, matrix, 0, matrixSize);
	}

	public void setNodePartialsForUpdate(int nodeIndex) {

		currentPartialsIndices[nodeIndex] = 1 - currentPartialsIndices[nodeIndex];
	}

	/**
	 * Sets the currently updating node partials for node nodeIndex. This may
	 * need to repeatedly copy the partials for the different category partitions
	 */
	public void setCurrentNodePartials(int nodeIndex, double[] partials) {

		double[] expandedPartials = partials;
		if (partials.length < partialsSize) {
			expandedPartials = new double[partials.length * matrixCount];
			int k = 0;
			for (int i = 0; i < matrixCount; i++) {
				System.arraycopy(partials, 0, expandedPartials, k, partials.length);
				k += partials.length;
			}
		}

		assert this.ptrPartials[0][nodeIndex] > DEBUG;
		assert partials != null;

		setNativeMemoryArray(expandedPartials, 0, this.ptrPartials[0][nodeIndex], 0, expandedPartials.length);
	}

	/**
	 * Calculates partial likelihoods at a node.
	 *
	 * @param nodeIndex1 the 'child 1' node
	 * @param nodeIndex2 the 'child 2' node
	 * @param nodeIndex3 the 'parent' node
	 */
	public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {

		if (ptrStates[nodeIndex1] != 0) {
			if (ptrStates[nodeIndex2] != 0) {
				calculateStatesStatesPruning(
						ptrStates[nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrStates[nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
//		        double[] tmp = new double[partialsSize];
//		        getNativeMemoryArray(ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3],0,tmp,0,partialsSize);
//		        System.err.println("NATIVE: STATE-STATE");
//		        System.err.println(new Vector(tmp));
			} else {
				calculateStatesPartialsPruning(
						ptrStates[nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrPartials[currentPartialsIndices[nodeIndex2]][nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
//		        double[] tmp = new double[partialsSize];
//		        getNativeMemoryArray(ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3],0,tmp,0,partialsSize);
//		        System.err.println("NATIVE: STATE-PARTIAL");
//		        System.err.println(new Vector(tmp));
			}
		} else {
			if (ptrStates[nodeIndex2] != 0) {
				calculateStatesPartialsPruning(
						ptrStates[nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex1]][nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
//			    double[] tmp = new double[partialsSize];
//		        getNativeMemoryArray(ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3],0,tmp,0,partialsSize);
//		        System.err.println("NATIVE: PARTIAL-STATE");
//		        System.err.println(new Vector(tmp));
			} else {
//			    System.err.println("Doing: PARTIAL-PARTIAL");
				calculatePartialsPartialsPruning(
						ptrPartials[currentPartialsIndices[nodeIndex1]][nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrPartials[currentPartialsIndices[nodeIndex2]][nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3]);
//			    double[] tmp = new double[partialsSize];
//		        getNativeMemoryArray(ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3],0,tmp,0,partialsSize);
//		        System.err.println("NATIVE: PARTIAL-PARTIAL");
//		        System.err.println(new Vector(tmp));
			}
		}

		if (useScaling) {
			scalePartials(nodeIndex3);
		}

//
//        int k =0;
//        for (int i = 0; i < patternCount; i++) {
//            double f = 0.0;
//
//            for (int j = 0; j < stateCount; j++) {
//                f += partials[currentPartialsIndices[nodeIndex3]][nodeIndex3][k];
//                k++;
//            }
//            if (f == 0.0) {
//                Logger.getLogger("error").severe("A partial likelihood (node index = " + nodeIndex3 + ", pattern = "+ i +") is zero for all states.");
//            }
//        }
	}

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(long ptrStates1, long ptrMatrices1,
	                                            long ptrStates2, long ptrMatrices2,
	                                            long ptrPartials3) {

		assert ptrStates1 > DEBUG;
		assert ptrMatrices1 > DEBUG;
		assert ptrStates2 > DEBUG;
		assert ptrMatrices2 > DEBUG;
		assert ptrPartials3 > DEBUG;

		nativeStatesStatesPruning(ptrStates1, ptrMatrices1, ptrStates2, ptrMatrices2, patternCount, matrixCount, ptrPartials3, stateCount);
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(long ptrStates1, long ptrMatrices1,
	                                              long ptrPartials2, long ptrMatrices2,
	                                              long ptrPartials3) {

		assert ptrStates1 > DEBUG;
		assert ptrMatrices1 > DEBUG;
		assert ptrPartials2 > DEBUG;
		assert ptrMatrices2 > DEBUG;
		assert ptrPartials3 > DEBUG;

		nativeStatesPartialsPruning(ptrStates1, ptrMatrices1, ptrPartials2, ptrMatrices2, patternCount, matrixCount, ptrPartials3, stateCount);
	}

	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(long ptrPartials1, long ptrMatrices1,
	                                                long ptrPartials2, long ptrMatrices2,
	                                                long ptrPartials3) {

		assert ptrPartials1 > DEBUG;
		assert ptrMatrices1 > DEBUG;
		assert ptrPartials2 > DEBUG;
		assert ptrMatrices2 > DEBUG;
		assert ptrPartials3 > DEBUG;

		nativePartialsPartialsPruning(ptrPartials1, ptrMatrices1, ptrPartials2, ptrMatrices2, patternCount, matrixCount, ptrPartials3, stateCount);
	}

	/**
	 * Calculates partial likelihoods at a node.
	 *
	 * @param nodeIndex1 the 'child 1' node
	 * @param nodeIndex2 the 'child 2' node
	 * @param nodeIndex3 the 'parent' node
	 * @param matrixMap  a map of which matrix to use for each pattern (can be null if integrating over categories)
	 */
	public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap) {
		if (ptrStates[nodeIndex1] != 0) {
			if (ptrStates[nodeIndex2] != 0) {
				calculateStatesStatesPruning(
						ptrStates[nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrStates[nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap);
			} else {
				calculateStatesPartialsPruning(
						ptrStates[nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrPartials[currentPartialsIndices[nodeIndex2]][nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap);

			}
		} else {
			if (ptrStates[nodeIndex2] != 0) {
				calculateStatesPartialsPruning(
						ptrStates[nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex1]][nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap);
			} else {
				calculatePartialsPartialsPruning(
						ptrPartials[currentPartialsIndices[nodeIndex1]][nodeIndex1], ptrMatrices[currentMatricesIndices[nodeIndex1]][nodeIndex1],
						ptrPartials[currentPartialsIndices[nodeIndex2]][nodeIndex2], ptrMatrices[currentMatricesIndices[nodeIndex2]][nodeIndex2],
						ptrPartials[currentPartialsIndices[nodeIndex3]][nodeIndex3], matrixMap);
			}
		}

		if (useScaling) {
			scalePartials(nodeIndex3);
		}
	}

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(long ptrStates1, long ptrMatrices1,
	                                            long ptrStates2, long ptrMatrices2,
	                                            long ptrPartials3, int[] matrixMap) {

		throw new RuntimeException("calculateStatesStatesPruning not implemented using matrixMap");
	}


	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(long ptrStates1, long ptrMatrices1,
	                                              long ptrPartials2, long ptrMatrices2,
	                                              long ptrPartials3, int[] matrixMap) {

		throw new RuntimeException("calculateStatesPartialsPruning not implemented using matrixMap");
	}


	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(long ptrPartials1, long ptrMatrices1,
	                                                long ptrPartials2, long ptrMatrices2,
	                                                long ptrPartials3, int[] matrixMap) {

		throw new RuntimeException("calculatePartialsPartialsPruning not implemented using matrixMap");

	}


	public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
		calculateIntegratePartials(ptrPartials[currentPartialsIndices[nodeIndex]][nodeIndex], proportions, outPartials);
	}

	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {

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


	/**
	 * Integrates partials across categories.
	 *
	 * @param ptrInPartials the partials at the node to be integrated
	 * @param proportions   the proportions of sites in each category
	 * @param outPartials   an array into which the integrated partials will go
	 */
//    protected abstract void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials);
	protected void calculateIntegratePartials(long ptrInPartials, double[] proportions, double[] outPartials) {

		assert ptrInPartials > DEBUG;

		nativeIntegratePartials(ptrInPartials, proportions, patternCount, matrixCount, outPartials, stateCount);
	}

	protected void scalePartials(int nodeIndex) {

		throw new UnsupportedOperationException("scaling partials is not yet available");

//        int u = 0;
//
//        for (int i = 0; i < patternCount; i++) {
//
//            double scaleFactor = 0.0;
//            int v = u;
//            for (int k = 0; k < matrixCount; k++) {
//                for (int j = 0; j < stateCount; j++) {
//                    if (partials[currentPartialsIndices[nodeIndex]][nodeIndex][v] > scaleFactor) {
//                        scaleFactor = partials[currentPartialsIndices[nodeIndex]][nodeIndex][v];
//                    }
//                    v++;
//                }
//                v += (patternCount - 1) * stateCount;
//            }
//
//            if (scaleFactor < scalingThreshold) {
//
//                v = u;
//                for (int k = 0; k < matrixCount; k++) {
//                    for (int j = 0; j < stateCount; j++) {
//                        partials[currentPartialsIndices[nodeIndex]][nodeIndex][v] /= scaleFactor;
//                        v++;
//                    }
//                    v += (patternCount - 1) * stateCount;
//                }
//                scalingFactors[currentPartialsIndices[nodeIndex]][nodeIndex][i] = Math.log(scaleFactor);
//
//            } else {
//                scalingFactors[currentPartialsIndices[nodeIndex]][nodeIndex][i] = 0.0;
//            }
//            u += stateCount;
//        }
	}

	/**
	 * This function returns the scaling factor for that pattern by summing over
	 * the log scalings used at each node. If scaling is off then this just returns
	 * a 0.
	 *
	 * @return the log scaling factor
	 */
	public double getLogScalingFactor(int pattern) {

		return 0;

//        double logScalingFactor = 0.0;
//        if (useScaling) {
//            for (int i = 0; i < nodeCount; i++) {
//                logScalingFactor += scalingFactors[currentPartialsIndices[i]][i][pattern];
//            }
//        }
//        return logScalingFactor;
	}

	/**
	 * Gets the partials for a particular node.
	 *
	 * @param nodeIndex   the node
	 * @param outPartials an array into which the partials will go
	 */
	public void getPartials(int nodeIndex, double[] outPartials) {

		assert ptrPartials[currentPartialsIndices[nodeIndex]][nodeIndex] > DEBUG;
		assert outPartials != null;

		getNativeMemoryArray(ptrPartials[currentPartialsIndices[nodeIndex]][nodeIndex], 0,
				outPartials, 0, partialsSize);
	}

	/**
	 * Store current state
	 */
	public void storeState() {

		System.arraycopy(currentMatricesIndices, 0, storedMatricesIndices, 0, nodeCount);
		System.arraycopy(currentPartialsIndices, 0, storedPartialsIndices, 0, nodeCount);
	}

	/**
	 * Restore the stored state
	 */
	public void restoreState() {
		// Rather than copying the stored stuff back, just swap the pointers...
		int[] tmp1 = currentMatricesIndices;
		currentMatricesIndices = storedMatricesIndices;
		storedMatricesIndices = tmp1;

		int[] tmp2 = currentPartialsIndices;
		currentPartialsIndices = storedPartialsIndices;
		storedPartialsIndices = tmp2;
	}

	/* Native memory handing functions */

	// Debug wrapper

	private long wrapperAllocateNativeMemoryArray(int length) {

		long rtnValue = allocateNativeMemoryArray(length);

		assert rtnValue > DEBUG;

		return rtnValue;
	}

	private void wrapperFreeNativeMemoryArray(long nativePtr) {
		assert nativePtr > DEBUG;

		freeNativeMemoryArray(nativePtr);
	}

	private native long allocateNativeMemoryArray(int length);

	private native void freeNativeMemoryArray(long nativePtr);

	private native void setNativeMemoryArray(double[] fromArray, int fromOffset,
	                                         long toNativePtr, int toOffset, int length);

	private native void setNativeMemoryArray(int[] fromArray, int fromOffset,
	                                         long toNativePtr, int toOffset, int length);

	private native void getNativeMemoryArray(long fromNativePtr, int fromOffset,
	                                         double[] toArray, int toOffset, int length);

	private native void getNativeMemoryArray(long fromNativePtr, int fromOffset,
	                                         int[] toArray, int toOffset, int length);

	/* Native peeling functions */


	public native void nativeIntegratePartials(long ptrPartials, double[] proportions,
	                                           int patternCount, int matrixCount,
	                                           double[] outPartials,
	                                           int stateCount);

	protected native void nativePartialsPartialsPruning(long ptrPartials1, long ptrMatrices1,
	                                                    long ptrPartials2, long ptrMatrices2,
	                                                    int patternCount, int matrixCount,
	                                                    long ptrPartials3,
	                                                    int stateCount);

	protected native void nativeStatesPartialsPruning(long ptrStates1, long ptrMatrices1,
	                                                  long ptrPartials2, long ptrMatrices2,
	                                                  int patternCount, int matrixCount,
	                                                  long ptrPartials3,
	                                                  int stateCount);

	protected native void nativeStatesStatesPruning(long ptrStates1, long ptrMatrices1,
	                                                long ptrStates2, long ptrMatrices2,
	                                                int patternCount, int matrixCount,
	                                                long ptrPartials3,
	                                                int stateCount);

	/* Library loading routines */

	public static boolean isAvailable() {
		return isNativeAvailable;
	}

	private static boolean isNativeAvailable = false;

	static {

		try {
			System.loadLibrary("NativeMemoryLikelihoodCore");
			isNativeAvailable = true;
		} catch (UnsatisfiedLinkError e) {
		}

	}

}