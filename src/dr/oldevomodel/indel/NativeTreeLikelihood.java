/*
 * NativeTreeLikelihood.java
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




/**
 * NativeTreeLikelihood - calls an implementation of HomologyRecursion.treeRecursion in C.
 *
 * @author Gerton Lunter
 *
 * 10/9/2003 
 *
 */

public class NativeTreeLikelihood {
	
  private static boolean isNativeAvailable = false;
  private int[][] iSequences;                       // rest are references to existing arrays
  private int iNN;
  private int iMUD;
  private long iStaticDataHandle;
    
    static {
	try {

	    System.loadLibrary("NativeTreeLikelihood");
	    System.err.println("Fast TKF91 homology integrator code found");
	    isNativeAvailable = true;

	} catch (UnsatisfiedLinkError e) {

	    System.err.println("Using Java TKF91 homology integrator code");

	}
    }
    

    /**
     * Constructor
     *
     */

    public NativeTreeLikelihood() { }


    /**
     * Hello, anybody?
     *
     */

    public boolean isAvailable() {
	return isNativeAvailable;
    }

    /**
     * Initialisation.  Returns with 'handle' (actually a pointer) to static data.
     *
     */

    private native long nativeInit(int[] iParent, double[] iEquil, double[] iTrans, 
				   double[] iH, double[] iN, double[] iE, double[] iB);

    /**
     * The actual recursion.  'iSD' is the handle returned by nativeInit.
     *
     */

    private native double nativeTreeRecursion(int[] iSignature, int[] iColumn, int iNumNucs, int iMUD, long iSD);

    /**
     * Relinquishes handle.
     *
     */

    private native void nativeDestruct(long iSD);


    /**
     * Initialiser
     *
     */

    public void init(int iNumNucs, int iMaxUnalignDimension, int[] iParent, double[] iEquil, double[][][] iTrans,
		     int[][] iSequences0,
		     double[] iN, double[] iH, double[] iE, double[] iB)
    {
	if (!isNativeAvailable)
	    return;

	// copy data into local variables
	iNN = iNumNucs;
        iMUD = iMaxUnalignDimension;
	iSequences = iSequences0;

	// translate iTrans[][][] array into flat array
	int iNumNodes = iTrans.length;
	double[] iFlatTrans = new double[ iNN*iNN*iNumNodes ];
	for (int i=0; i<iNumNodes; i++) {
	    for (int j=0; j<iNN; j++) {
		for (int k=0; k<iNN; k++) {
		    iFlatTrans[ k + iNN*( j + iNN*i ) ] = iTrans[i][j][k];
		}
	    }
	}

	// store data someplace where the C code has fast access
	iStaticDataHandle = nativeInit(iParent, iEquil, iFlatTrans, iH, iN, iE, iB);
    }


    /**
     * Make sure that memory is released when this object is garbage collected
     *
     */

    protected void finalize() throws Throwable {

	if (iStaticDataHandle != 0) {
	    nativeDestruct( iStaticDataHandle );
	}
	super.finalize();
    }


    /**
     * Actual front-end for HomologyRecursion
     *
     */

    public double treeRecursion(IntMathVec iSignature, IntMathVec iPos)
    {

	int[] iColumn = new int[ iSignature.iV.length ];
	for (int i=0; i<iSignature.iV.length; i++) {
	    if (iSignature.iV[i] != 0)
		iColumn[i] = iSequences[i][iPos.iV[i]];
	}

	return nativeTreeRecursion(iSignature.iV,
				   iColumn,
				   iNN,
				   iMUD,
				   iStaticDataHandle);
    }

}	

	

