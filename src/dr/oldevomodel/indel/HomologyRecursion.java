/*
 * HomologyRecursion.java
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

package dr.oldevomodel.indel;


import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.math.BFloat;

import java.util.HashMap;

public class HomologyRecursion {
    
    /** Array of 0-1 vectors, encoding the alignment as columns */
    IntMathVec[] iAlignment;
    
    /** Array of arrays, encoding the sequences (in rows).  No gaps.  Range [0,iNumNucs).  */
    int[][] iSequences;

    /** 
     *   Array encoding the tree.  
     *   iParent[iChild] is the parent of iChild, except -1 = root's parent.
     *   Parent must be > child, and leaf indices must correspond to indices into alignment vectors
     */
    private int[] iParent;

    /** Length of incoming branches to each tree node.  Root is assumed to be infinity */
    private double[] iTau;

    /** Parameters of the TKF91 model */
    double iLambda;
    double iMu;

    /** Number of nucleotides */
    int iNumNucs;

    /** Transition probability matrices, iTrans[node][from][to], for each tree node except root */
    double[][][] iTrans;

    /** Equilibrium probabilities for nucleotides. */
    double[] iEquil;

    /** Intermediate results for TKF91 model, initialised in recursion(), used in treeRecursion() */
    private double iH[], iN[], iB[], iE[];

    /** TKF91-related result for initialising the DP table */
    private double iInitial;
    
    /** Class that implements a fast version of treeLikelihood */
    private NativeTreeLikelihood iNativeMethod;

    /** Labels used in local array iState[] */
    static private final int eFree = 0, ePossible = 1, eEdgeUsed = 2, eUsed = 3;

    /** Max dimension of unalignable region */
    static private final int cMaxUnalignDimension = 10;
    
    /** Minimum edge length */
    static private final double MIN_EDGE_LENGTH = 1e-3;

    /** Error count */
    static int sBigUnalignableRegion = 0;

    String PrintDouble(double[] d) {
	String str = "";
        for(double aD : d) {
            str = str + aD + " ";
        }
	return str;
    }

    void checkConsistency() {
		// Checks whether variables are init'ed properly
		int iCols = iAlignment.length;
		int iRows = iAlignment[0].iV.length;
		String iErr = "";
		IntMathVec iPos = new IntMathVec( iRows );
		for (int i=0; i<iCols; i++) {
		    for (int j=0; j<iRows; j++) {
			if (iAlignment[i].iV[j] != 0 && iAlignment[i].iV[j] != 1) {
			    iErr = "Non-0/1 emissions.";
			}
		    }
		    iPos.add( iAlignment[i] );
		}
		for (int j=0; j<iRows; j++) {
		    if (iSequences[j].length != iPos.iV[j]) {
			iErr = "Bad sequences length";
		    }
		    for (int i=0; i<iPos.iV[j]; i++) {
			if ((iSequences[j][i]<0)||(iSequences[j][i]>=iNumNucs)) {
			    iErr = "Nucleotide codes in iSequences not in range, " + iSequences[j][i] + "not in [0,"+(iNumNucs-1)+"]";
			}
		    }
		}
		int[] iChildren = new int[ iParent.length ];
		for (int i=0; i<iParent.length-1; i++)
		    iChildren[ iParent[i] ]++;
		for (int i=0; i<iRows; i++)
		    if (iChildren[i] != 0) {
			iErr = "Num tips does not correspond to num sequences, or bad tree";
		    }
		for (int i=iRows; i<iParent.length; i++)
		    if (iChildren[i] != 2) {
			iErr = "Bad tree - not binary?  Or too many tips?";
		    }
		if (iParent[iParent.length-1] != -1) {
		    iErr = "Bad tree - root not in final position, or not labeled -1";
		}
		for (int i=0; i<iParent.length-1; i++) {
			if (iTau[i] < MIN_EDGE_LENGTH) {
				iErr = "";//Bad tree - edge lengths below threshold (edge " + i + ", len " + iTau[i] + ")";
				iTau[i] = MIN_EDGE_LENGTH;
			}
		}
		
		//for (int i=0; i < iAlignment.length; i++) {
		//	System.out.println(iAlignment[i]);
		//}

		if ( !iErr.equals("") )
		    System.out.println(iErr);
	}    



    public HomologyRecursion() {
	// Nothing?  Initialisation is done manually through init() below.
    }
    
    public void init(Tree tree, Alignment alignment, SubstitutionModel substModel, 
		     double mutationRate, double lengthDistr, double deathRate) {
   
    	// initialize the iParent and iTau arrays based on the given tree.
	initTree(tree, mutationRate);
	
	int[] treeIndex = new int[tree.getTaxonCount()];
	for (int i =0; i < treeIndex.length; i++) {
	    treeIndex[i] = tree.getTaxonIndex(alignment.getTaxonId(i));
	    //System.out.println("alignment[" + i + "] = tree[" + treeIndex[i] + "]");
	}
	  		
    	// initialize the iAlignment array from the given alignment.
    	initAlignment(alignment, treeIndex);
    	
    	// initialize the iSequences array from the given alignment.
    	initSequences(alignment, treeIndex);
    	
    	// initialize the iTrans array from the substitution model -- must be called after populating tree!
    	initSubstitutionModel(substModel);
    	
    	// iLambda, iMu
    	iLambda = deathRate*lengthDistr;
    	iMu = deathRate;
    	
    	// iNumNucs - alphabet size 
	DataType dataType = substModel.getDataType();
    	iNumNucs = dataType.getStateCount();
	
	// Initialise TKF91 coefficients in iB, iH, iN, iE, and iInitial
	initTKF91();

	// Check
	checkConsistency();

	// Initialise native method
	iNativeMethod = new NativeTreeLikelihood();
	iNativeMethod.init(iNumNucs, cMaxUnalignDimension, iParent, iEquil, iTrans, iSequences, iN, iH, iE, iB);

    }
    
    private void initTree(Tree tree, double mutationRate) {
	iParent = new int[tree.getNodeCount()];
	iTau = new double[tree.getNodeCount()-1];
	populate(tree, tree.getRoot(), new int[] {tree.getExternalNodeCount()}, mutationRate);
	iParent[tree.getNodeCount()-1] = -1;
	
    } 
    
   /** 
    * initialize the iTrans array from the substitution model -- must be called after populating tree!
    */
    private void initSubstitutionModel(SubstitutionModel model) {
    	
    	DataType dataType = model.getDataType();
    	int stateCount = dataType.getStateCount();
    	
    	iTrans = new double[iTau.length][stateCount][stateCount];
    	
    	double[] transProb = new double[stateCount*stateCount];
    	int count; 
    	for (int i =0; i < iTau.length; i++) {
    		model.getTransitionProbabilities(iTau[i], transProb);
    		count = 0;
    		for (int j = 0; j < stateCount; j++) {
    			for (int k = 0; k < stateCount; k++) {
    				iTrans[i][j][k] = transProb[count];
    				count += 1;
    			}
    		}
    	}
    	
    	// initialize equlibrium distribution
    	iEquil = new double[stateCount];
    	for (int k = 0; k < stateCount; k++) {
			iEquil[k] = model.getFrequencyModel().getFrequency(k);
		}
    }
    
    /**
     * Initializes the iAlignment array from the given alignment.
     */
    private void initAlignment(Alignment alignment, int[] treeIndex) {
    	
    	int numSeqs = alignment.getSequenceCount();
    	int numSites = alignment.getSiteCount();
    	DataType dataType = alignment.getDataType();
    	int numStates = dataType.getStateCount();
    	
    	iAlignment = new IntMathVec[numSites];
    	
    	int[] column = new int[numSeqs];
    	for (int i =0; i < numSites; i++) {
    		for (int j = 0; j < numSeqs; j++) {
   				column[treeIndex[j]] = ((alignment.getState(j, i) >= numStates) ? 0 : 1);	
    		}
    		iAlignment[i] = new IntMathVec(column);
    	}
    }
    
     /**
     * Initializes the iSequence array from the given alignment.
     */
    private void initSequences(Alignment alignment, int[] treeIndex) {
    	
    	int numSeqs = alignment.getSequenceCount();
    	DataType dataType = alignment.getDataType();
    	int numStates = dataType.getStateCount();
    	
    	iSequences = new int[numSeqs][];
    	
    	for (int i = 0; i < numSeqs; i++) {
    	
    		int seqLength = 0;
    		for (int j =0; j < alignment.getSiteCount(); j++) {
				int state = alignment.getState(i, j);
				if (state>=0 && state<numStates) { seqLength += 1;}
				
			}
			iSequences[treeIndex[i]] = new int[seqLength];
			int count = 0;
			for (int j =0; j < alignment.getSiteCount(); j++) {
				int state = alignment.getState(i, j);
				if (state>=0 && state<numStates) {
					iSequences[treeIndex[i]][count] = state;
					count += 1;
				}
			}
    	}
    }
    
    /**
     * Populates the iParent and iTau arrays.
     * @return the node number of the given node.
     */
    private int populate(Tree tree, NodeRef node, int[] current, double mutationRate) {
    	
    	int nodeNumber = node.getNumber();
    	
    	// if its an external node just return the number
    	if (tree.isExternal(node)) {
    		iTau[nodeNumber] = 
    			(tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node))  * mutationRate; 
    		return nodeNumber;
    	}
    	
    	// if internal node, first let your children be assigned numbers
    	int[] childNumbers = new int[tree.getChildCount(node)];
    	for (int i = 0; i < tree.getChildCount(node); i++) {
    		childNumbers[i] = populate(tree, tree.getChild(node, i), current, mutationRate);
    	}
    	
    	// now, pick the next available number
    	nodeNumber = current[0];
    	// if you are not the root, then record the branch length above you.
    	if (!tree.isRoot(node)) { 
    		//iTau[nodeNumber] = tree.getBranchLength(node) * mutationRate; 
    		iTau[nodeNumber] = 
    			(tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node))  * mutationRate; 
    	}
    	// increment the next available number
    	current[0] += 1;
    	
    	// now that you have your number, populate the iParent entries of your children.
    	for (int i = 0; i < tree.getChildCount(node); i++) {
    		iParent[childNumbers[i]] = nodeNumber;
    	}
    	
    	// finally return your number so your parent can do the same.
    	return nodeNumber;
    }


    /**
     * Populates the iB, iE, iH and iN arrays, and iInitial value
     * @return void
     */
  public void initTKF91() {
    
    int iNumNodes = iParent.length;
    double iBeta[] = new double[iNumNodes];
    iB = new double[iNumNodes];
    iE = new double[iNumNodes];
    iH = new double[iNumNodes];
    iN = new double[iNumNodes];
    iInitial = 1.0;
    for (int i=0; i<iNumNodes; i++) {
      if (i==iNumNodes-1) {
	// root
	iBeta[ i ] = 1.0/iMu;
	iH[i] = 0.0;
      } else {
	// internal node or tip
	iBeta[ i ] = Math.exp((iLambda-iMu)*iTau[i]);
	iBeta[ i ] = (1.0-iBeta[i])/(iMu-iLambda*iBeta[i]);
	iH[i] = Math.exp(-iMu*iTau[i])*(1.0-iLambda*iBeta[i]);
      }
      iB[i] = iLambda*iBeta[i];
      iE[i] = iMu*iBeta[i];
      iN[i] = (1.0-iMu*iBeta[i])*(1.0-iB[i]) - iH[i];
      iInitial *= (1.0 - iB[i]);
    }
  }
	

    /** 
     * Calculates 'conditional factor' associated to emission encoded by iSignature, at position iPos.
     * Nulls in iSignature represent no-emissions, others indicate homology family (1,2,...,cMaxUnalignDimension, inclusive)
     */
    private double treeRecursion(IntMathVec iSignature, IntMathVec iPos) {

	// See if we can use C code
	if (iNativeMethod.isAvailable()) {
	    return iNativeMethod.treeRecursion(iSignature, iPos);
	}

	int iLeaves = iSignature.iV.length;          // Dimension of alignment columns, i.e. number of leaves
	int iNumNodes = iParent.length;           // Number of internal nodes
	int[] iHom = new int[ iNumNodes ];        // Homology for every node; 0 if node need not be homologous to an emitted nucleotide
	int[] iHomNum = new int[ iNumNodes ];     // Number of homologous emissions accounted for by homologous nucleotide @ this node
	int[] iHomMultiplicity = new int[ cMaxUnalignDimension+1 ];   	   // Number of emissions for each class of homologous nucleotides
	int[] iChild1 = new int[ iNumNodes ], iChild2 = new int[ iNumNodes ];     	   // Left and right children
	double[][] iFh = new double[ iNumNodes ][],  iFi = new double[ iNumNodes ][];      // Fhomolgous and Finhomologous arrays

	// Initialise stuff
	for (int i=0; i<iNumNodes; i++) {
	    iFh[i] = new double[ iNumNucs ];
	    iFi[i] = new double[ iNumNucs + 1];   // Extra position for 'gap' entry
	}
	for (int i=0; i<iLeaves; i++) {
	    iHomMultiplicity[ iSignature.iV[ i ] ]++;
	    iHom[ i ] = iSignature.iV[ i ];
	    if (iSignature.iV[i] == 0)
		iHomNum[ i ] = 0;
	    else
		iHomNum[ i ] = 1;
	}

	// Loop over all nodes except root, and find out which nodes need carry nucleotides of what 
	// homology class.  Also fill iChild* arrays, which point to the two children of every parent.
	boolean iClashingHomology = false;
	for (int i=0; i<iNumNodes-1; i++) {
	    if ((iHomNum[i] == iHomMultiplicity[ iHom[i] ]) || (iHom[i] == 0)) {
			// This node is the MRCA of the homologous nucleotides iHom[i], or a gap - do nothing
	    } else {
			if (iHom[ iParent[ i ]] == 0 || iHom[ iParent[ i ]] == iHom[i]) {
			    // This node is not yet MRCA, so node above must be homologous
			    iHom[iParent[i]] = iHom[i];
			    // If iHom[i]==0, iHomNum[i]==0.
			    iHomNum[iParent[i]] += iHomNum[i];
			} else {
			    // Clashing homology; signal
			    iClashingHomology = true;
			}
	    }
	    if (iChild1[ iParent[ i ] ] == 0)
		iChild1[ iParent[ i ] ] = i;
	    else
		iChild2[ iParent[ i ] ] = i;
	}

	// Bail out - cheaper than do this implicitly in the recursion below
	if (iClashingHomology) {
	    return 0.0;
	}

	// Start recursion.  First initialise the leaves
	for (int i=0; i<iLeaves; i++) {
	    if (iSignature.iV[i] == 0) {
		// gap
		iFi[i][iNumNucs] = 1.0;
	    } else {
		// nucleotide
		iFh[i][ iSequences[i][iPos.iV[i]] ] = 1.0;
	    }
	    //System.out.println("Fi" + i +" = " + PrintDouble(iFi[i]));
	    //System.out.println("Fh" + i +" = " + PrintDouble(iFh[i]));
	}

	// Now do the recursion, bottom-up, on all internal nodes
	for (int i=iLeaves; i<iNumNodes; i++) {
	    // Find out whether:
	    // 1- One homology family spanning tree intersects both child edges, i.e.
	    //    a 'homologous' nucleotide should travel down both edges,
	    // 2- One homology family Spanning tree intersects one of the child edges, i.e.
	    //    'homologous' nucleotide must travel down a specific edge
	    // 3- No homology family's spanning tree intersect either child edge, i.e.
	    //    a 'homologous' nucleotide may travel down either edge and pop out at a leaf,
	    //    or may travel down both edges but has to die in at least one subtree, or
	    //    an 'inhomologous' nucleotide may do whatever it likes here.
	    if ((iHom[i] != 0) && (iHom[i] == iHom[iChild1[i]]) && (iHom[i] == iHom[iChild2[i]])) {
			// case 1

			for (int j=0; j<iNumNucs; j++) {
			    double iL = 0.0;
			    double iR = 0.0;
			    for (int k=0; k<iNumNucs; k++) {
				iL += iFh[iChild1[i]][k] * iH[iChild1[i]] * iTrans[iChild1[i]][j][k];
				iR += iFh[iChild2[i]][k] * iH[iChild2[i]] * iTrans[iChild2[i]][j][k];
			    }
			    iFh[i][j] = iL*iR;
			}
			// Others: 0.0
	    } else if (iHom[i] != 0) {
			// case 2.  Figure out which is the homologous child, and which the inhomologous one.

			int iChildH, iChildI;
			if (iHom[i] == iHom[iChild1[i]]) {
			    iChildH = iChild1[i];
			    iChildI = iChild2[i];
			} else {
			    iChildH = iChild2[i];
			    iChildI = iChild1[i];
			}
	
			for (int j=0; j<iNumNucs; j++) {
			    double iL = 0.0;
			    double iR = iE[iChildI] * iFi[iChildI][iNumNucs];
			    for (int k=0; k<iNumNucs; k++) {

					iL += iFh[iChildH][k] * iH[iChildH] * iTrans[iChildH][j][k];
					iR += 
				    	(iFh[iChildI][k] + iFi[iChildI][k]) * (iN[iChildI] - iE[iChildI]*iB[iChildI]) * iEquil[k] +
				    	iFi[iChildI][k] * iH[iChildI] * iTrans[iChildI][j][k];

			    	}

			    iFh[i][j] = iL*iR;

			}
			// Others: 0.0
	
	    } else {
			// case 3
			int iC1 = iChild1[i];
			int iC2 = iChild2[i];
			//System.out.println(iC1);
			//System.out.println(iC2);
			for (int j=0; j<iNumNucs; j++) {
			    double iL1 = 0.0, iL2 = 0.0;
			    double iR1 = iE[iC1] * iFi[iC1][iNumNucs];
			    double iR2 = iE[iC2] * iFi[iC2][iNumNucs];
			    //System.out.println(iR1);
			    //System.out.println(iR2);
			    for (int k=0; k<iNumNucs; k++) {
			    	iL1 += iFh[iC1][k] * iH[iC1] * iTrans[iC1][j][k];
					iL2 += iFh[iC2][k] * iH[iC2] * iTrans[iC2][j][k];
					iR1 += 
					    (iFh[iC1][k] + iFi[iC1][k]) * (iN[iC1] - iE[iC1]*iB[iC1]) * iEquil[k] +
					    iFi[iC1][k] * iH[iC1] * iTrans[iC1][j][k];
					iR2 += 
					    (iFh[iC2][k] + iFi[iC2][k]) * (iN[iC2] - iE[iC2]*iB[iC2]) * iEquil[k] +
					    iFi[iC2][k] * iH[iC2] * iTrans[iC2][j][k];
			    }
			    iFh[i][j] = iL1*iR2 + iL2*iR1;  // homology pops out below iC1 + homology pops out below iC2
			    iFi[i][j] = iR1*iR2;            // no homology with j below i.
			}
	
			double iL = iFi[iC1][iNumNucs];
			double iR = iFi[iC2][iNumNucs];
			for (int j=0; j<iNumNucs; j++) {
		    
			    iL -= iB[iC1] * (iFi[iC1][j] + iFh[iC1][j]) * iEquil[j];
				iR -= iB[iC2] * (iFi[iC2][j] + iFh[iC2][j]) * iEquil[j];
	
			}
			iFi[i][iNumNucs] = iL*iR;
	    }  

	}   //  recursion over the internal nodes

	// Now calculate the final result

	int iRoot = iNumNodes-1;
	double iResult = iFi[iRoot][iNumNucs];
	for (int i=0; i<iNumNucs; i++)
	    iResult -= (iFi[iRoot][i] + iFh[iRoot][i]) * iB[iRoot] * iEquil[i];

	return iResult;

    }  // treeRecursion


	/** 
	 * @return the logLikelihood of the alignment given tree, substitution model, et cetera
	 */
	 
    public double recursion() {

	// Initialise

	int iLen = iAlignment.length;
	int iLeaves = iAlignment[0].iV.length;
	int iFirstNotUsed = 0;                         // First not-'used' alignment vector (for efficiency) 
	int iState[] = new int[ iLen ];                // Helper array, to traverse the region in the DP table corresp. to the alignment
	IntMathVec iPos = new IntMathVec( iLeaves );   // Current position; sum of all used vectors
        HashMap<IntMathVec, BFloat> iTable = new HashMap<IntMathVec, BFloat>();

	// Calculate correction factor for null emissions ("wing folding", or linear equation solving.)

	double iNullEmissionFac = treeRecursion( iPos, iPos );
	//System.out.println("Null emisison: " + iNullEmissionFac);

	// Enter first probability into DP table
	iTable.put( iPos, new BFloat(iInitial / iNullEmissionFac) );

	// Array of possible vector indices, used in inner loop 
	int[] iPossibles = new int[cMaxUnalignDimension];

	do {

	    // Find all possible vectors from current position, iPos
	    IntMathVec iMask = new IntMathVec( iLeaves );
	    int iPtr;
	    int iNumPossible = 0;
	    for (iPtr = iFirstNotUsed; iMask.zeroEntry() && iPtr<iLen; iPtr++) {
			if (iState[ iPtr ] != eUsed) {
		    	if (iMask.innerProduct( iAlignment[iPtr] ) == 0) {
					iState[ iPtr ] = ePossible;
					//System.out.println("Accepting column " + iAlignment[iPtr] + " as number " + iNumPossible);
					if (iNumPossible == cMaxUnalignDimension) {
			    		// This gets too hairy - bail out.
			    		sBigUnalignableRegion++;
			    		System.err.println("We bailed out cause it was hairy: iNumPossible=" + iNumPossible);
			    		return Double.NEGATIVE_INFINITY;
					}
					iPossibles[iNumPossible++] = iPtr;
		    	}	
		    	iMask.add( iAlignment[iPtr] );
			}
	    }

	    // Loop over all combinations of possible vectors, which define edges from 
	    // iPos to another possible position, by ordinary binary counting.

	    IntMathVec iNewPos = new IntMathVec( iPos );
	    IntMathVec iSignature = new IntMathVec( iPos.iV.length );
	    int iPosPtr;
	    boolean iUnusedPos;
	    boolean iFoundNonZero;

	    do {

			// Find next combination
			iFoundNonZero = false;
			for (iPosPtr = iNumPossible - 1; iPosPtr >= 0; --iPosPtr) {
			    int iCurPtr = iPossibles[ iPosPtr ];
			    if (iState[ iCurPtr ] == ePossible) {
					iState[ iCurPtr ] = eEdgeUsed;
					iNewPos.add( iAlignment[ iCurPtr ] );
					// Compute signature vector
					iSignature.addMultiple( iAlignment[ iCurPtr ], iPosPtr+1 );
					// Signal: non-zero combination found, and stop
					iFoundNonZero = true;
					iPosPtr = 0;
			    } else {
					// It was eEdgeUsed (i.e., digit == 1), so reset digit and continue
					iState[ iCurPtr ] = ePossible;
					iNewPos.subtract( iAlignment[ iCurPtr ] );
					iSignature.addMultiple( iAlignment[ iCurPtr ], -iPosPtr-1 );
			    }
			}

			if (iFoundNonZero) {
			    //System.out.print("Reading from pos " + iPos);
			    BFloat iLeft = (BFloat)(iTable.get( iPos )).clone();
			    //System.out.print(" left=" + iLeft);
			    BFloat iRightObj = iTable.get( iNewPos );
			    BFloat iRight;
			    if (iRightObj == null) {
					iUnusedPos = true;
					iRight = new BFloat(0);
			    } else {
					iRight = iRightObj;
					iUnusedPos = false;
			    }
			    double iTransFac = (-treeRecursion( iSignature, iPos )) / iNullEmissionFac;
			    //System.out.print(" fac=" + iTransFac);
			    iLeft.multiply( iTransFac );
			    iRight.add( iLeft ); 
			    //System.out.print(" sig=" + iSignature);
			    // And store
			    //System.out.println(" Storing pos " + iNewPos+ " val=" + iRight);
			    // If we are storing a value at a previously unused position, make sure we use a fresh key object
			    if (iUnusedPos) {
					iTable.put( iNewPos.clone(), iRight );
			    } else {
					iTable.put( iNewPos, iRight );
			    }
			}
	
		} while (iFoundNonZero);
	
		// Now find next entry in DP table.  Use farthest unused vector
		--iPtr;
		while (iPtr >= 0 && iState[iPtr] != ePossible) {
			// Undo any possible used vector that we encounter
			if (iState[iPtr] == eUsed) {
		    	iPos.subtract( iAlignment[iPtr] );
		    	iState[iPtr] = eFree;
			}
			--iPtr;
	    }
	    
	    if (iPtr == -1) {
			// No more unused vectors, so we also fell through the edge loop above,
			// hence iNewPos contains the final position
			//System.out.println("Returning " + (BFloat)iTable.get( iNewPos ) + " = " + ((BFloat)iTable.get( iNewPos ) ).log());
			return (iTable.get( iNewPos )).log();
	    }

	    // Now use this farthest-out possible vector
	    iState[iPtr] = eUsed;
	    iPos.add( iAlignment[iPtr] );
	    if (iPtr <= iFirstNotUsed)
			iFirstNotUsed++;
	    
	} while (true);

    } // recursion
}

