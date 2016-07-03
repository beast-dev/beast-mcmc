/*
 * IstvansProposal.java
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

import java.util.HashMap;
import java.util.Random;


public class IstvansProposal {
    
    static double cGap = -5;
    // Gap penalty
    static int cGapsymbol = -1;
    // Is -1 representing gap?
    static double cBasis = 1.5;

    static final int cMaxlength = 1000;
    
    static final int cMaxUnalignDimension = 17;
    
    static final int eFree = 0;
    static final int ePossible = 1;
    static final int eEdgeUsed = 2;
    static final int eUsed = 3;
    
    // Static count of path count bound overruns
    static int sBigUnalignableRegion = 0;
    private final Random r = new Random();
    // dynamic programming table
    // does this need initialization before each proposal?
    private final double[][] iDP = new double[cMaxlength][cMaxlength];

    public void setGapSymbol(int gapSymbol) {
		cGapsymbol = gapSymbol;
    }
    
    /**
     * @returns the combinatorial factor associated to the alignment, i.e.
     * the number of valid paths in the DP table.
     */
    int countPaths(int[][] iInputAlignment, int iStartCol, int iEndCol) {
	
		// Initialise
		
		int iLen = iEndCol - iStartCol + 1;
		int iLeaves = iInputAlignment.length;
		int iFirstNotUsed = 0;                         // First not-'used' alignment vector (for efficiency) 
		int iState[] = new int[ iLen ];                // Helper array, to traverse the region in the DP table corresp. to the alignment
		IntMathVec iPos = new IntMathVec( iLeaves );   // Current position; sum of all used vectors
        HashMap<IntMathVec, Integer> iTable = new HashMap<IntMathVec, Integer>();
		
		IntMathVec[] iAlignment = new IntMathVec[iLen];
		for (int i=iStartCol; i<=iEndCol; i++) {
		    iAlignment[i-iStartCol] = new IntMathVec(iLeaves);
		    for (int j=0; j<iLeaves; j++) {
			iAlignment[i-iStartCol].iV[j] = ((iInputAlignment[j][i] == cGapsymbol) ? 0 : 1);
		    }
		}
		
		// Enter first count into DP table
		iTable.put(iPos, 1);
		
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
				if (iNumPossible == cMaxUnalignDimension) {
				    // This gets too hairy - bail out.
				    sBigUnalignableRegion++;
				    return -1;
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
			    int iLeft = iTable.get(iPos);
			    //System.out.print(" left=" + iLeft);
			    int iRight;
			    Object iRightObj = iTable.get( iNewPos );
			    if (iRightObj == null) {
				iRight = 0;
				iUnusedPos = true;
			    } else {
				iRight = (Integer) iRightObj;
				iUnusedPos = false;
			    }
			    //System.out.print(" right=" + iRight);
			    iRight += iLeft;
			    //System.out.print(" sig=" + iSignature);
			    // And store
			    //System.out.println(" Storing pos " + iNewPos);
			    // If we are storing a value at a previously unused position, make sure we use a fresh key object
			    if (iUnusedPos) {
				iTable.put(iNewPos.clone(), iRight);
			    } else {
				iTable.put( iNewPos, iRight);
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
			return iTable.get(iNewPos);
		    }
		    
		    
		    // Now use this farthest-out possible vector
		    iState[iPtr] = eUsed;
		    iPos.add( iAlignment[iPtr] );
		    if (iPtr <= iFirstNotUsed)
			iFirstNotUsed++;
		    
		} while (true);
		
    } // recursion
    
    
    
    /**
     * @returns the log of the hastings ratio
     * I cannot guarantee it works, however :) Istvan 12/06/2003
     * Version 2.0 It returns with a statistically correct Hastings ratio
     */
    public double propose(int[][] iAlignment, double[][][] iProbs, double[] iBasefreqs, int[] iParent, double[] iTau, int[][] returnedAlignment, double iP, double exponent, double gapPenalty) {
	
		int           iL1;            
		// length of the cut out alignment
		int           iPos;
		// position of the cut out alignment (first index)
		int           iL2;            
		// length of the obtained alignment
		int[][][]     iProfile = new int[iParent.length][iParent.length][cMaxlength];      
		// profiles obtained in the proposal
		int[][]       iTempProfile = new int[iParent.length][cMaxlength];
		// temporary profile since at the stochastic traceback we get the profile in reverse order
		int[]         iXP = new int[iParent.length], iYP = new int[iParent.length];       
		// length of a profile (X) and number of sequences involved into the profile (Y)
		double        iTempdouble1 = 0.0, iTempdouble2 = 0.0, iTempdouble3 = 0.0, iTempdouble4 = 0.0;
		// temporal variables in the dp recursion and stochastic traceback
		double        iProposalLoglikelihood = 0.0;
		// loglikelihood of the proposal, what else?
		int           iTempint1;
		// temporal variables for counting gaps, iTempint1 is also used at child-sorting
		double[][]    iScoring = new double[iProbs[0].length][iProbs[0][0].length];
		//  scoring table for the dynamic programming
		int[]         iChild1 = new int[iParent.length], iChild2 = new int[iParent.length];
		// these arrays make the work on the tree easier
		int[]         iSorter = new int[iParent.length];
		// this array helps at sorting children so the profile of root will be the correct multiple alignment
		boolean       iTempBoole;
		// temporary boolean variable for checking all-gap columns
		int[][]       iProfileNumber = new int[iParent.length][iAlignment.length];
		
		
		//Testing preconditions
		if (iP <= 0.0 || iP >= 1.0) throw new IllegalArgumentException("iP must be in range (0,1)");
		if (exponent < 1.0) throw new IllegalArgumentException("exponent must be in range [1,+infinity)");
		if (gapPenalty >= 0.0) throw new IllegalArgumentException("gapPenalty must be in range (-infinity,0)");
		
		cGap = gapPenalty;
		cBasis = exponent;
		
		
		
		/* obtaining iChild1 and iChild2 */
		
		for(int i = 0; i < iParent.length - 1; i++) {
		    if(iChild1[iParent[i]] == 0) {
			iChild1[iParent[i]] = i;
		    } else {
			iChild2[iParent[i]] = i;
		    }
		}
		
		// now sorting children...
		
		for(int i = 0; i < iAlignment.length; i++)
		    iSorter[i]=i;
		// iSorter for leaves
		for(int i = iAlignment.length; i < iParent.length; i++) {
		    if(iSorter[iChild1[i]] > iSorter[iChild2[i]])
			iSorter[i] = iSorter[iChild2[i]];
		    else
			iSorter[i] = iSorter[iChild1[i]];
		    // iSorter for internal nodes, it gets the lower value
		    
		    if(iSorter[iChild1[i]] > iSorter[iChild2[i]]) {
			// iChild1[i] has the smaller iSorter value
			// this guarantees that the root's profile contains the sequences in the proper order
			iTempint1 = iChild1[i];
			iChild1[i] = iChild2[i];
			iChild2[i] = iTempint1;
		    }
		}
		/* cutting a part of the alignment */
		
		// obtaining the length of the cut out part (iL1) and the first index (iPos)
		for(iL1=1; iL1 < iAlignment[0].length && Math.random() > iP; iL1++);
		
		
		//iL1 = 1;
		iPos = r.nextInt(iAlignment[0].length - iL1 + 1);
		
		//System.out.println("Window size = " + iL1);
		// collecting the cut out sequences from the alignment
		// they will be the first profiles
		for(int i = 0; i < iAlignment.length; i++) {
		    iYP[i] = 1;
		    for(int j = iPos; j <iPos + iL1; j++)
			if(iAlignment[i][j] != cGapsymbol) {
			    iProfile[i][0][iXP[i]] = iAlignment[i][j];
			    iXP[i]++;
			}
		    // iProfileNumber
		    iProfileNumber[i][0] = i;
		}
		
		//Printing out the cut out sequences
		//	System.out.println("Cut out sequences:");
		//for(int i = 0; i < iAlignment.length; i++){
		//   System.out.print(i+". ");
		//  for(int j = 0; j < iXP[i]; j++){
		//System.out.print(iProfile[i][0][j]);
		//  }
		//  System.out.println("");
		//}
		
		
		/* dynamic programming algorithm for the proposal */

		for(int k = iAlignment.length; k < iParent.length; k++) {
		    
		    // obtaining the scoring matrix
		    for(int i = 0; i < iProbs[0].length; i++) {
				for(int j = 0; j < iProbs[0][0].length; j++){
				    iScoring[i][j] = 0;
				    for(int l = 0; l < iProbs[0].length; l++)
						iScoring[i][j] += iProbs[iChild1[k]][i][l] * iProbs[iChild2[k]][l][j];
				}
			}
		    for(int i = 0; i < iProbs[0].length; i++) {
				for(int j = 0; j < iProbs[0][0].length; j++) {
				    //  System.out.println("iScoring["+i+"]["+j+"] = "+iScoring[i][j]);
				    iScoring[i][j] = Math.log(iScoring[i][j] / iBasefreqs[j]);
				     	//System.out.println("iScoring["+i+"]["+j+"] = "+iScoring[i][j]);
				}
		    }
		    
		    // initialising the table
		    iDP[0][0]=0.0;
		    for(int i = 1; i <= iXP[iChild1[k]]; i++){
			iDP[i][0] = iDP[i - 1][0] + cGap;
		    }
		    for(int j = 1; j <= iXP[iChild2[k]]; j++) {
			iDP[0][j] = iDP[0][j - 1] + cGap;
		    }
		    
		    // main dynamic programming
		    for(int i = 1; i <= iXP[iChild1[k]]; i++)
			for(int j = 1; j <= iXP[iChild2[k]]; j++){
			    iTempdouble1 = iDP[i - 1][j] + cGap;
			    
			    iTempdouble2 = iDP[i][j - 1] + cGap;
			    
			    {
				int iCounter = 0;
			    iTempdouble3 = 0.0;
			    for(int l = 0; l < iYP[iChild1[k]]; l++)
				for(int m = 0; m < iYP[iChild2[k]]; m++) {
				    if(iProfile[iChild1[k]][l][i - 1] != cGapsymbol && iProfile[iChild2[k]][m][j - 1] != cGapsymbol){
					iTempdouble3 += iScoring[iProfile[iChild1[k]][l][i - 1]][iProfile[iChild2[k]][m][j - 1]];
					iCounter++;
				    }
				}
			    iTempdouble3 = iDP[i - 1][j - 1] + iTempdouble3/iCounter;
			    }
			
			    
			    if(iTempdouble1 > iTempdouble2)
				iDP[i][j] = iTempdouble1;
			    else
				iDP[i][j] = iTempdouble2;
			    if(iTempdouble3 > iDP[i][j])
				iDP[i][j] = iTempdouble3;
			}
		    
		    // Stochastic traceback, generating the new profile
		    int iI = iXP[iChild1[k]];
		    int jJ = iXP[iChild2[k]];
		    while(iI > 0 || jJ > 0) {
				if(iI > 0){
				    iTempdouble1 = iDP[iI - 1][jJ] + cGap;
				}
				
				
				if(jJ > 0) {
				    iTempdouble2 = iDP[iI][jJ - 1] + cGap;
				}
				
				if(iI > 0 && jJ > 0) {
				    int iCounter = 0;
				    iTempdouble3 = 0.0;
				    for(int l = 0; l < iYP[iChild1[k]]; l++)
					for(int m = 0; m < iYP[iChild2[k]]; m++) {
					    if(iProfile[iChild1[k]][l][iI - 1] != cGapsymbol && iProfile[iChild2[k]][m][jJ - 1] != cGapsymbol){
						iTempdouble3 += iScoring[iProfile[iChild1[k]][l][iI - 1]][iProfile[iChild2[k]][m][jJ - 1]];
						iCounter++;
					    }
					}
				    iTempdouble3 = iDP[iI - 1][jJ - 1] + iTempdouble3/iCounter;
				}
				
				//System.out.println("iTempdouble1 = " + iTempdouble1);
				//System.out.println("iTempdouble2 = " + iTempdouble2);
				//System.out.println("iTempdouble3 = " + iTempdouble3);
				
				
				// normalising
				iTempdouble4 = 0.0;
				if(iI > 0) {
				    iTempdouble4 += Math.exp(iTempdouble1 * Math.log(cBasis)); 
				}
				if(jJ > 0) {
				    iTempdouble4 += Math.exp(iTempdouble2 * Math.log(cBasis));
				}
				if(iI > 0 && jJ > 0) {
				    iTempdouble4 += Math.exp(iTempdouble3 * Math.log(cBasis));
				}
				
				//System.out.println("iTempdouble4 = " + iTempdouble4);
				
				
				if(iI > 0)
				    iTempdouble1 = Math.exp(iTempdouble1 * Math.log(cBasis))/iTempdouble4;
				else
				    iTempdouble1 = 0.0;
				
				if(jJ > 0)
				    iTempdouble2 = Math.exp(iTempdouble2 * Math.log(cBasis))/iTempdouble4;
				else
				    iTempdouble2 = 0.0;
				
				if(iI > 0 && jJ > 0)
				    iTempdouble3 = Math.exp(iTempdouble3 * Math.log(cBasis))/iTempdouble4;
				else
				    iTempdouble3 = 0.0;

				//System.out.println("Sum: " + (iTempdouble1+iTempdouble2+iTempdouble3));
				
				
				// stochastic step
				iTempdouble4 = Math.random();
				if(iTempdouble4 < iTempdouble1) {
					// the probability is in the denominator, that's why we substract
				    iProposalLoglikelihood -= Math.log(iTempdouble1);
				    for(int l = 0; l < iYP[iChild1[k]]; l++) {
					iTempProfile[l][iXP[k]] = iProfile[iChild1[k]][l][iI - 1];
				    }
				    for(int l = iYP[iChild1[k]]; l < iYP[iChild1[k]] + iYP[iChild2[k]]; l++) {
					iTempProfile[l][iXP[k]] = cGapsymbol;
				    }         
				    iXP[k]++;
				    iI--;
				    //	    if(iI == -1) {  // it must almost never happen, however, if it does happen in a magical way...
				    //	iXP[k]--;
				    //	iI++;
				    //}
				} 
				else if(iTempdouble4 < iTempdouble1 + iTempdouble2) {
					// the probability is in the denominator, that's why we substract
				    iProposalLoglikelihood -= Math.log(iTempdouble2);
				    for(int l = 0; l < iYP[iChild1[k]]; l++) {
					iTempProfile[l][iXP[k]] = cGapsymbol;
				    }    
				    for(int l = iYP[iChild1[k]]; l < iYP[iChild1[k]] + iYP[iChild2[k]]; l++) {
					iTempProfile[l][iXP[k]] = iProfile[iChild2[k]][l - iYP[iChild1[k]]][jJ - 1];
				    }
				    iXP[k]++;
				    jJ--;
				    //if(jJ == -1) {  // it must almost never happen, however, if it does happen in a magical way...
				    //	iXP[k]--;
				    //	jJ++;
				    //}
				}
				else {
				    if(iTempdouble3 <= 0) {
						System.out.println("How the fuck could it happen?");
				    }
					// the probability is in the denominator, that's why we substract
				    iProposalLoglikelihood -= Math.log(iTempdouble3);
				    for(int l = 0; l < iYP[iChild1[k]]; l++)
					iTempProfile[l][iXP[k]] = iProfile[iChild1[k]][l][iI - 1];
				    for(int l = iYP[iChild1[k]]; l < iYP[iChild1[k]] + iYP[iChild2[k]]; l++)
					iTempProfile[l][iXP[k]] = iProfile[iChild2[k]][l - iYP[iChild1[k]]][jJ - 1];
				    iXP[k]++;
				    jJ--;
				    iI--;
				    //   if(jJ == -1 || iI == -1) {  // it must almost never happen, however, if it does happen in a magical way...
				    //	iXP[k]--;
				    //	iI++;
				    //	jJ++;
				    //}
				    
				}
			
			
			
		    }// finished the traceback, iI, jJ
		    
		    
		    //System.out.println("iProposalLoglikelihood:"+iProposalLoglikelihood);
		    
		    iYP[k] = iYP[iChild1[k]] + iYP[iChild2[k]];
		    
		    //Making new iProfileNumbers
		    for(int i = 0; i < iYP[iChild1[k]]; i++)
				iProfileNumber[k][i] = iProfileNumber[iChild1[k]][i];
		    for(int i = iYP[iChild1[k]]; i < iYP[iChild1[k]] + iYP[iChild2[k]]; i++)
				iProfileNumber[k][i] = iProfileNumber[iChild2[k]][i - iYP[iChild1[k]]];
		    
		    // now reversing the profile, loading from the iTempProfile into iProfile
		    
		    for(int i = 0; i < iXP[k]; i++) {
				for(int j = 0; j < iYP[k]; j++) {
					iProfile[k][j][i] = iTempProfile[j][iXP[k] - i - 1];
				}
			}
		    
		} // k, index for the internal nodes, iProfile[root] is the obtained multiple alignment

		iL2 = iXP[iParent.length - 1];
		
		// putting the new sub-alignment into the new proposal
		
		//System.out.println("iL1 = " + iL1);
		//System.out.println("iL2 = " + iL2);
		//System.out.println("iPos = " + iPos);
		//System.out.println("iParent.length = " + iParent.length);
		//System.out.println("iAlignment.length = " + iAlignment.length);
		//System.out.println("returnedAlignment.length = " + returnedAlignment.length);
		//System.out.println("iAlignment[0].length = " + iAlignment[0].length);
		
		
		
		for (int i = 0; i < iAlignment.length; i++) {
		    returnedAlignment[i] = new int[iAlignment[0].length + iL2 - iL1];
		    //System.out.println("returnedAlignment[i].length = " + returnedAlignment[i].length);
		    
		}
		
		for(int i = 0; i < iPos; i++) {
		    for(int j = 0; j < iAlignment.length; j++) {
			returnedAlignment[j][i] = iAlignment[j][i];
		    }
		}
		for(int i = iPos; i < iPos + iL2; i++) {
		    for(int j = 0; j < iAlignment.length; j++) {
			//System.out.println("i = " + i + ", j = " + j);
		     	returnedAlignment[iProfileNumber[iParent.length - 1][j]][i] = iProfile[iParent.length - 1][j][i - iPos];
		    }
		}
		
		
		//Printing the new alignment
		//System.out.println("Proposed alignment");
		//	for(int i = 0; i < iAlignment.length; i++){
		//  System.out.print(i+". ");
		//   for(int j = 0; j < iL2; j++){
		//System.out.print(iProfile[iParent.length-1][i][j]);
		//  }
		//  System.out.println("");
		//}
		
		
		
		
		for(int i = iPos + iL1; i < iAlignment[0].length; i++) {
		    for(int j = 0; j < iAlignment.length; j++) {
			returnedAlignment[j][i + iL2 -iL1] = iAlignment[j][i];
		    }
		}
		
		
		// I involve the combinatorical factor into the Hastings ratio
		// countPaths(int[][] iInputAlignment, int iStartCol, int iEndCol)
		
		
		int iCF = countPaths(returnedAlignment,iPos,iPos+iL2-1);
		if(iCF == -1){
		    System.out.println("****   HEY ALEXEI! IT'S A FUCKING FUCKED ALIGNMENT! ****");
		    iProposalLoglikelihood -= Math.exp(100.0);
		    iCF = 1;
		}
		
		
		iProposalLoglikelihood -= Math.log((double)iCF);
		
		
		iCF = countPaths(iAlignment,iPos,iPos+iL1-1);
		if(iCF == -1){
		    System.out.println("****   HEY ALEXEI! IT'S A FUCKING FUCKED ALIGNMENT! ****");
		    iProposalLoglikelihood -= Math.exp(100.0);
		    iCF = 1;
		}
		iProposalLoglikelihood += Math.log((double)iCF);
		
		
		/*
		  Now I'm going to calculate the Hastings ratio.
		  To do that, I must make the profile of each node
		  I copy the multiple alignment into profiles, and then delete the all-gaps columns
		  The traceback is easy: based on checking all-gap subcolumns
		*/
		
		
		// collecting the sequences from the alignment -- now without cut out!
		// they will be the first profiles
		for(int i = 0; i < iAlignment.length; i++) {
		    iYP[i] = 1;
		    for(int j = iPos; j <iPos + iL1; j++)
			iProfile[i][0][j - iPos] = iAlignment[i][j];
		}
		
		
		// profiles at internal nodes are obtained by merging profiles of children
		for(int i = iAlignment.length; i < iParent.length; i++) {
		    iYP[i] = iYP[iChild1[i]] + iYP[iChild2[i]];
		    for(int j = 0; j < iL1; j++) {
			for(int k = 0; k < iYP[iChild1[i]]; k++)
			    iProfile[i][k][j] = iProfile[iChild1[i]][k][j];
			for(int k = iYP[iChild1[i]]; k < iYP[iChild1[i]] + iYP[iChild2[i]]; k++)
			    iProfile[i][k][j] = iProfile[iChild2[i]][k - iYP[iChild1[i]]][j];
		    }
		}
		
		
		// and now -- cutting out all-gap columns
		for(int i = 0; i < iParent.length; i++){
		    // first i load it into the temporary profile and then back...
		    for(int j = 0; j < iL1; j++)
			for(int k = 0; k < iYP[i]; k++)
			    iTempProfile[k][j] = iProfile[i][k][j];
		    
		    
		    iXP[i] = 0;
		    for(int j = 0; j < iL1; j++){
				// checking for all-gap;
				boolean notAllGaps = false;
				for(int k = 0; k < iYP[i]; k++){
					if (iTempProfile[k][j] != cGapsymbol) { 
						notAllGaps = true;
						break;
					}
				}

				if(notAllGaps) {
				    // if not all-gap then 
				    for(int k = 0; k < iYP[i]; k++)
					iProfile[i][k][iXP[i]] =  iTempProfile[k][j];
				    iXP[i]++;
				}
		    }
	       
		}
	
	
	// profiles are ready, so now the DP
	for(int k = iAlignment.length; k < iParent.length; k++) {
	    
	    
	    // obtaining the scoring matrix
	    for(int i = 0; i < iProbs[0].length; i++)
		for(int j = 0; j < iProbs[0][0].length; j++){
		    iScoring[i][j] = 0.0;
		    for(int l = 0; l < iProbs[0].length; l++)
			iScoring[i][j] += iProbs[iChild1[k]][i][l] * iProbs[iChild2[k]][l][j];
		}
	    for(int i = 0; i < iProbs[0].length; i++)
		for(int j = 0; j < iProbs[0][0].length; j++)
		    iScoring[i][j] = Math.log(iScoring[i][j] / iBasefreqs[j]);
	    
	    
	    // initialising the table
	    iDP[0][0]=0.0;
	    for(int i = 1; i <= iXP[iChild1[k]]; i++){
		iDP[i][0] = iDP[i - 1][0] + cGap;
	    }
	    for(int j = 1; j <= iXP[iChild2[k]]; j++) {
		iTempdouble2 = 0.0;
		for(int m = 0; m < iYP[iChild2[k]]; m++)
		    if(iProfile[iChild2[k]][m][j - 1] != cGapsymbol)
			iTempdouble2 += cGap;
		iDP[0][j] = iDP[0][j - 1] + cGap;
	    }
	    
	    
	    // main dynamic programming
	    for(int i = 1; i <= iXP[iChild1[k]]; i++)
		for(int j = 1; j <= iXP[iChild2[k]]; j++){
		    iTempdouble1 = iDP[i - 1][j] + cGap;
		    
		    
		    iTempdouble2 = 0.0;
		    for(int m = 0; m < iYP[iChild2[k]]; m++)
			if(iProfile[iChild2[k]][m][j - 1] != cGapsymbol)
			    iTempdouble2 += cGap;
		    iTempdouble2 = iDP[i][j - 1] + cGap;
		    
		    {
			int iCounter = 0;
		    iTempdouble3 = 0.0;
		    for(int l = 0; l < iYP[iChild1[k]]; l++)
			for(int m = 0; m < iYP[iChild2[k]]; m++) {
			    if(iProfile[iChild1[k]][l][i - 1] != cGapsymbol && iProfile[iChild2[k]][m][j - 1] != cGapsymbol){
				iTempdouble3 += iScoring[iProfile[iChild1[k]][l][i - 1]][iProfile[iChild2[k]][m][j - 1]];
				iCounter++;
			    }
			}
		    iTempdouble3 = iDP[i - 1][j - 1] + iTempdouble3/iCounter;
		    }
		    
		    if(iTempdouble1 > iTempdouble2)
			iDP[i][j] = iTempdouble1;
		    else
			iDP[i][j] = iTempdouble2;
		    if(iTempdouble3 > iDP[i][j])
			iDP[i][j] = iTempdouble3;
		}
	    
	    
	    // traceback, now not stochastic
	    int i = iXP[iChild1[k]];
	    int j = iXP[iChild2[k]];
	    int n = iXP[k];
	    while(n > 0) {
		
		
		if(i > 0){
		    //System.out.println("i = " + i + "; j = " + j);
		    
		    iTempdouble1 = iDP[i - 1][j] + cGap;
		}
		
		
		
		if(j > 0) {
		    //System.out.println("i = " + i + "; j = " + j);
		     
		    iTempdouble2 = iDP[i][j - 1] + cGap;
		}
		
		
		if(i > 0 && j > 0) {
		    int iCounter = 0;
		    iTempdouble3 = 0.0;
		    for(int l = 0; l < iYP[iChild1[k]]; l++)
			for(int m = 0; m < iYP[iChild2[k]]; m++) {
			    if(iProfile[iChild1[k]][l][i - 1] != cGapsymbol && iProfile[iChild2[k]][m][j - 1] != cGapsymbol){
				iTempdouble3 += iScoring[iProfile[iChild1[k]][l][i - 1]][iProfile[iChild2[k]][m][j - 1]];
				iCounter++;
			    }
			}
		    iTempdouble3 = iDP[i - 1][j - 1] + iTempdouble3/iCounter;
		}
		
		
		// normalising
		iTempdouble4 = 0.0;
		if(i > 0)
		    iTempdouble4 += Math.exp(iTempdouble1 * Math.log(cBasis)); 
		if(j > 0)
		    iTempdouble4 += Math.exp(iTempdouble2 * Math.log(cBasis));
		if(i > 0 && j > 0)
		    iTempdouble4 += Math.exp(iTempdouble3 * Math.log(cBasis));
		
		
		if(i > 0)
		    iTempdouble1 = Math.exp(iTempdouble1 * Math.log(cBasis))/iTempdouble4;
		else
		    iTempdouble1 = 0.0;
		
		
		if(j > 0)
		    iTempdouble2 = Math.exp(iTempdouble2 * Math.log(cBasis))/iTempdouble4;
		else
		    iTempdouble2 = 0.0;
		
		
		if(i > 0 && j > 0)
		    iTempdouble3 = Math.exp(iTempdouble3 * Math.log(cBasis))/iTempdouble4;
		else
		    iTempdouble3 = 0.0;
		
		//System.out.println("Sum in the second DP: " + (iTempdouble1+iTempdouble2+iTempdouble3));
		
		
		// finding the proper step (all-gap subcolumns tell us the way
		iTempBoole = true;
		for(int l = 0; l < iYP[iChild1[k]] && iTempBoole; l++)
		    iTempBoole = (iProfile[k][l][n-1] == cGapsymbol);
		
		
		if(iTempBoole) {
		    // all the first part is gap, so we decrease j
		    j--;
		    iProposalLoglikelihood += Math.log(iTempdouble2);
		    n--;
		}
		else {
		    iTempBoole = true;
		    for(int l = iYP[iChild1[k]]; l < iYP[iChild1[k]] + iYP[iChild2[k]] && iTempBoole; l++)
			iTempBoole = (iProfile[k][l][n-1] == cGapsymbol);
		    if(iTempBoole) {
			// all the second part is gap, so we decrease i
			i--;
			iProposalLoglikelihood += Math.log(iTempdouble1);
			n--;
		    }
		    else {
			// no all-gap subcolumn, decreasing both i and j;
			i--;
			j--;
			iProposalLoglikelihood += Math.log(iTempdouble3);
			n--;
		    }
		}
		
		
	    } // end of traceback
	    
	    //System.out.println("iProposalLogLikelihood2=" + iProposalLoglikelihood);
	} // end of k
	
	//System.out.println("iL1=" + iL1);
	//System.out.println("iL2=" + iL2);
	return iProposalLoglikelihood + (iL1 - iL2)*Math.log(1.0-iP);
    }
}