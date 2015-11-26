/*
 * DiffusionLikelihood.java
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

package dr.evolution.continuous;

import dr.math.*;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

/**
 *
 *
 */
public class DiffusionLikelihood implements UnivariateFunction, MultivariateFunction {

    Tree tree;
    String traitName;

    public DiffusionLikelihood(Tree tree, String traitName) {
        this.tree = tree;
        this.traitName = traitName;
    }

    public double optimize(double[] p) {
        //

        if (p.length > 1) {
            p[0] = 1.0;
            p[1] = 0.01;

            MultivariateMinimum optimizer = new ConjugateDirectionSearch();
            //MultivariateMinimum optimizer = new DifferentialEvolution(2,5);
            optimizer.optimize(this,p,1e-6, 1e-6);
            return -evaluate(p);
        } else if (p.length == 1) {
            UnivariateMinimum optimizer = new UnivariateMinimum();
            p[0] = optimizer.optimize(this,1e-6);
            return -evaluate(p[0]);
        } else throw new IllegalArgumentException("");
    }

    /**
     * @return the log likelihood of going from start to stop in the given time
     */
    public double getLogLikelihood(double distance, double time, double D, double bias) {

        // expected variance of distances of given time
        double Dtime = D * time;

        double d = distance - (bias*time);


        //System.out.println("distance=" + unbiasedDistance + " time=" + time);

        // the log likelihood of travelling distance d, in time t given diffusion rate D
        return -(Math.log(Math.sqrt(Dtime*2*Math.PI))) - ((d*d)/(2*Dtime));
    }


    public double evaluate(double[] argument) {

        double D = argument[0];
        double bias = argument[1];

        double logLkl = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {

            NodeRef child = tree.getNode(i);

            if (!tree.isRoot(child)) {

                NodeRef parent = tree.getParent(child);

                Contrastable parentValue = (Contrastable)tree.getNodeAttribute(parent,traitName);
                Contrastable childValue = (Contrastable)tree.getNodeAttribute(child,traitName);

                double distance = parentValue.getDifference(childValue);
                //System.out.println(distance);

                logLkl += getLogLikelihood(distance,tree.getBranchLength(child),D, bias);

            }
        }

        return -logLkl;
    }




    public int getNumArguments() {
        return 2;
    }

    public double getLowerBound(int n) {
        if (n == 0) return 1e-12;
        return -Double.MAX_VALUE;
    }

    public double getUpperBound(int n) {
        return 2000;
    }

    public double evaluate(double argument) {
        return evaluate(new double[] {argument, 0.0});
    }

    public double getLowerBound() {
        return getLowerBound(0);
    }

    public double getUpperBound() {
        return getUpperBound(0);
    }
}
