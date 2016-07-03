/*
 * MDSTest.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.multidimensionalscaling;

import dr.inference.multidimensionalscaling.NativeMDSSingleton;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class MDSTest extends MathTestCase {

    private static NativeMDSSingleton mds = loadLibrary();

    private static NativeMDSSingleton loadLibrary() {
        try {
            return NativeMDSSingleton.loadLibrary();
        } catch (UnsatisfiedLinkError error) {
            System.err.println("Unable to load MDS library; no trying tests");
            return null;
        }
    }

    public void testInitialization() {
        if (mds != null) {
            int i = mds.initialize(2, 100, 0);
            assertEquals(i, 0);
            i = mds.initialize(2, 100, 0);
            assertEquals(i, 1);
        } else {
            System.out.println("testInitialization skipped");
        }

    }

    public void testMakeDirty() {
        if (mds != null) {
            mds.makeDirty(0);
        } else {
            System.out.println("testMakeDirty skipped");
        }

    }


}
