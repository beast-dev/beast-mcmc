/*
 * ContourWithR.java
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

package dr.geo.contouring;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;

/**
 * @author Marc Suchard
 */
public class ContourWithR implements ContourMaker {


    public ContourWithR(final double[] xValues, final double[] yValues) {
        this(xValues, yValues, 50);
    }

    public ContourWithR(final double[] xValues, final double[] yValues, int N) {
        this.xValues = xValues;
        this.yValues = yValues;
        this.N = N;
    }

    private final double[] xValues;
    private final double[] yValues;

    public ContourPath[] getContourPaths(double hpdValue) {
        makeContour(xValues, yValues, hpdValue, N);
        if (contourList == null)
            return null;

        ContourPath[] contourPaths = new ContourPath[getNumberContours()];

        for(int i=0; i<getNumberContours(); i++) {
            double[][] cont = getContour(i);
            contourPaths[i] = new ContourPath(null,1,cont[0],cont[1]);
        }
        return contourPaths;
    }


    public void makeContour(double[] xValues, double[] yValues, double hpd) {
        makeContour(xValues, yValues, hpd, N);
    }

    public void makeContour(double[] xValues, double[] yValues, double hpd, int N) {


        REXP x = rEngine.eval("makeContour(" +
                makeRString(xValues) + "," +
                makeRString(yValues) + "," +
                hpd + "," +
                N + ")");
        contourList = x.asVector();
    }

    public int getNumberContours() {
        if (contourList != null)
            return contourList.size();
        return 0;
    }




    public double[][] getContour(int whichContour) {

        if (contourList != null) {
            double[][] result = new double[2][];
            RVector oneContour = contourList.at(whichContour).asVector();
            result[0] = oneContour.at(1).asDoubleArray();
            result[1] = oneContour.at(2).asDoubleArray();
            return result;
        }
        return null;
    }


    private static Rengine rEngine = null;
    private int N;


    private RVector contourList = null;

    private static final String[] rArgs = {"--no-save", "--max-vsize=1G"};

    private static final String[] rBootCommands = {
            "library(MASS)",
            "makeContour = function(var1, var2, prob=0.95, n=50, h=c(1,1)) {" +
                    "post1 = kde2d(var1, var2, n = n, h=h); " +    // This had h=h in argument
                    "dx = diff(post1$x[1:2]); " +
                    "dy = diff(post1$y[1:2]); " +
                    "sz = sort(post1$z); " +
                    "c1 = cumsum(sz) * dx * dy; " +
                    "levels = sapply(prob, function(x) { approx(c1, sz, xout = 1 - x)$y }); " +
                    "line = contourLines(post1$x, post1$y, post1$z, level = levels); " +
                    "return(line) }"
    };

    private String makeRString(double[] values) {
        StringBuffer sb = new StringBuffer("c(");
        sb.append(values[0]);
        for (int i = 1; i < values.length; i++) {
            sb.append(",");
            sb.append(values[i]);
        }
        sb.append(")");
        return sb.toString();
    }


    static public boolean processWithR = false;

    static {
        try {
            System.loadLibrary("jri");
            processWithR = true;
            System.err.println("JRI loaded. Will process using R contouring.");

//        if (!Rengine.versionCheck()) {
//            throw new RuntimeException("JRI library version mismatch");
//        }

        rEngine = new Rengine(rArgs, false, null);

        if (!rEngine.waitForR()) {
            throw new RuntimeException("Cannot load R");
        }

        for (String command : rBootCommands) {
            rEngine.eval(command);
        }
        } catch (UnsatisfiedLinkError e) {
            System.err.println("JRI not available. Will process using Java contouring.");
        }
    }
}
