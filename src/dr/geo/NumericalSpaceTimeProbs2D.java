/*
 * NumericalSpaceTimeProbs2D.java
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

package dr.geo;

import dr.math.distributions.MultivariateNormalDistribution;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * @author Alexei Drummond
 */
public class NumericalSpaceTimeProbs2D {

    final int latticeWidth;
    final int latticeHeight;
    final int tsteps;
    final int subtsteps;
    final double minx, miny, dx, dy, dt;
    final MultivariateNormalDistribution D;
    final SpaceTimeRejector rejector;

    int[][][][][] counts;
    int[][][] normalization;
    int[][][] maxCount;

    public NumericalSpaceTimeProbs2D(
            int latticeWidth,
            int latticeHeight,
            int tsteps,
            int subtsteps,
            double dt,
            Rectangle2D bounds,
            MultivariateNormalDistribution distribution,
            SpaceTimeRejector rejector) {

        this.latticeWidth = latticeWidth;
        this.latticeHeight = latticeHeight;
        this.tsteps = tsteps;
        this.subtsteps = subtsteps;
        this.D = distribution;
        this.rejector = rejector;
        this.dt = dt;

        minx = bounds.getMinX();
        miny = bounds.getMinY();

        dx = (bounds.getMaxX() - minx) / latticeWidth;
        dy = (bounds.getMaxY() - miny) / latticeHeight;

        counts = new int[latticeWidth][latticeHeight][latticeWidth][latticeHeight][tsteps];
        normalization = new int[latticeWidth][latticeHeight][tsteps];
        maxCount = new int[latticeWidth][latticeHeight][tsteps];
    }

    public void populate(Point2D start, int paths, boolean includeSubpaths) {

        populate(x(start.getX()), y(start.getY()), paths, includeSubpaths);
    }

    public int populateAbsorbing(Point2D start, int paths) {
        return populateAbsorbing(x(start.getX()), y(start.getY()), paths);
    }

    /**
     * @param i
     * @param j
     * @param paths
     * @return the number of successfully simulated paths
     */
    public int populateAbsorbing(int i, int j, int paths) {

        double subdt = dt / (double) subtsteps;

        double[] next = new double[2];
        double[] start = new double[2];
        int[] pathx = new int[tsteps];
        int[] pathy = new int[tsteps];
        int successes = 0;

        for (int reps = 0; reps < paths; reps += 1) {

            double time = 0.0;
            start[0] = (i + Math.random()) * dx + minx;
            start[1] = (j + Math.random()) * dy + miny;
            while (rejector.reject(0, start)) {
                start[0] = (i + Math.random()) * dx + minx;
                start[1] = (j + Math.random()) * dy + miny;
            }

            boolean reject = false;
            for (int t = 0; t < tsteps && !reject; t++) {

                for (int s = 0; s < subtsteps && !reject; s++) {

                    D.nextScaledMultivariateNormal(start, subdt, next);
                    time += subdt;
                    reject = rejector.reject(time, next);

                    if (!reject) {
                        start[0] = next[0];
                        start[1] = next[1];
                    }
                }

                if (!reject) {
                    pathx[t] = x(next[0]);
                    pathy[t] = y(next[1]);

                    increment(i, j, pathx[t], pathy[t], t);
                }
            }
            if (!reject) successes += 1;

            if (reps % 10000 == 0) {
                System.out.print(".");
                System.out.flush();
            }
        }
        System.out.println();
        return successes;
    }


    public void populate(int i, int j, int paths, boolean includeSubpaths) {

        double subdt = dt / (double) subtsteps;

        double[] next = new double[2];
        double[] start = new double[2];
        int[] pathx = new int[tsteps];
        int[] pathy = new int[tsteps];

        for (int reps = 0; reps < paths; reps += 1) {

            double time = 0.0;
            start[0] = (i + Math.random()) * dx + minx;
            start[1] = (j + Math.random()) * dy + miny;
            while (rejector.reject(0, start)) {
                start[0] = (i + Math.random()) * dx + minx;
                start[1] = (j + Math.random()) * dy + miny;
            }

            for (int t = 0; t < tsteps; t++) {

                for (int s = 0; s < subtsteps; s++) {

                    do {
                        D.nextScaledMultivariateNormal(start, subdt, next);
                        time += subdt;
                    } while (rejector.reject(time, next));

                    start[0] = next[0];
                    start[1] = next[1];
                }

                pathx[t] = x(next[0]);
                pathy[t] = y(next[1]);

                increment(i, j, pathx[t], pathy[t], t);
            }

            if (includeSubpaths) {
                for (int t = 0; t < tsteps; t++) {
                    for (int s = t + 1; s < tsteps; s++) {
                        increment(pathx[t], pathy[t], pathx[s], pathy[s], s - t - 1);
                    }
                }
            }

            if (reps % 1000 == 0) {
                System.out.print(".");
                System.out.flush();
            }
        }
    }

    private void increment(int i, int j, int k, int l, int t) {
        counts[i][j][k][l][t] += 1;
        normalization[i][j][t] += 1;

        if (counts[i][j][k][l][t] > maxCount[i][j][t]) {
            maxCount[i][j][t] = counts[i][j][k][l][t];
        }
    }

    public void populate(int paths) {

        System.out.println("Populating numerical transition probabilities");

        for (int i = 0; i < latticeWidth; i++) {

            for (int j = 0; j < latticeHeight; j++) {

                populate(i, j, paths, true);
            }
            System.out.print(".");
            System.out.flush();
        }

        System.out.println(latticeWidth * latticeHeight * paths + " new paths computed.");
    }

    public final int x(double x) {
        return (int) ((x - minx) / dx);
    }

    public final int y(double y) {
        return (int) ((y - miny) / dy);
    }

    public final int t(double time) {
        return (int) (time / dt);
    }

    public double getProb(Point2D start, Point2D end, double time) {

        int i = x(start.getX());
        int j = x(start.getY());
        int k = x(end.getX());
        int l = x(end.getY());

        if (time > tsteps * dt) {
            System.err.println("Time = " + time + ", max time estimated is " + tsteps * dt);
            return (double) counts[i][j][k][l][tsteps - 1] / normalization[i][j][tsteps - 1];
        } else {

            //time interpolation

            int t = t(time);
            double tlow = t * dt;
            double thigh = tlow + dt;

            double weightlow = (thigh - time) / dt;

            return weightlow * p(i, j, k, l, t) + (1.0 - weightlow) * p(i, j, k, l, t + 1);
        }
    }

    public double p(int i, int j, int k, int l, int t) {
        return (double) counts[i][j][k][l][t] / (double) normalization[i][j][t];

    }

    public double r(int i, int j, int k, int l, int t) {
        return (double) counts[i][j][k][l][t] / (double) maxCount[i][j][t];
    }

    public void writeToFile(String s) throws FileNotFoundException {

        PrintWriter writer = new PrintWriter(s);

        writer.write("xsteps=" + latticeWidth + "\n");
        writer.write("ysteps=" + latticeHeight + "\n");
        writer.write("tsteps=" + tsteps + "\n");
        writer.write("dx=" + dx + "\n");
        writer.write("dy=" + dy + "\n");
        writer.write("dt=" + dt + "\n");
        writer.write("minx=" + minx + "\n");
        writer.write("miny=" + miny + "\n");
        writer.write("D=" + matrixString());

        for (int i = 0; i < latticeWidth; i++) {
            for (int j = 0; j < latticeHeight; j++) {
                for (int k = 0; k < latticeWidth; k++) {
                    for (int l = 0; l < latticeHeight; l++) {
                        for (int t = 0; t < tsteps; t++) {

                            writer.write(i + "\t" + j + "\t" + k + "\t" + l + "\t" + t + "\t" + counts[i][j][k][l][t] + "\n");
                        }
                    }
                }
            }
        }
        writer.close();
    }

    private String matrixString() {

        double[][] m = D.getScaleMatrix();
        StringBuilder builder = new StringBuilder();
        builder.append("[");
        for (int i = 0; i < m.length; i++) {
            builder.append("[");
            builder.append(m[i][0]);
            for (int j = 1; j < m[i].length; j++) {
                builder.append("," + m[i][0]);
            }
            builder.append("]");

        }
        builder.append("]");
        return builder.toString();
    }

    public static void main(String[] args) throws FileNotFoundException {

        Rectangle2D bounds = new Rectangle2D.Double(0, 0, 1, 1);

        MultivariateNormalDistribution D =
                new MultivariateNormalDistribution(new double[]{0.0}, new double[][]{{1, 0}, {0, 1}});

        NumericalSpaceTimeProbs2D nstp2D =
                new NumericalSpaceTimeProbs2D(50, 50, 50, 1, 0.02, bounds, D, SpaceTimeRejector.Utils.createSimpleBounds2D(bounds));


        long startTime = System.currentTimeMillis();
        nstp2D.populate(0, 0, 1000, true);
        long stopTime = System.currentTimeMillis();
        System.out.println("Time taken = " + (stopTime - startTime) / 1000 + " seconds");


//        System.out.println("Writing to file...");
//        nstp2D.writeToFile("unitSquareDiffusion.txt");

        for (int i = 0; i < 10; i++) {

            Point2D start = new Point2D.Double(Math.random(), Math.random());
            Point2D end = new Point2D.Double(Math.random(), Math.random());
            double time = Math.random();

            double p = nstp2D.getProb(start, end, time);

            System.out.println("Pr(" + end.getX() + ", " + end.getY() + " | " + start.getX() + ", " + start.getY() + ", t=" + time + ") = " + p);
        }

    }

}
