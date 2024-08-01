/*
 * TimeVaryingBranchRateModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.math.distributions.GaussianProcessBasisApproximation;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.RandomFieldDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static dr.math.distributions.GaussianProcessBasisApproximation.getSpectralDensityEigenValue;
import static dr.math.distributions.GaussianProcessBasisApproximation.getSpectralDensity;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */
// TO DO: Pass GaussianProcessBasisApproximation as the parameter to remove code duplication
// TO DO: Make IntegratedSquaredGPApproximation a Random Field
public class IntegratedSquaredGPApproximation {


    private final double[] coefficients;
    private final double origin;
    private final double boundary;
    private final double marginalVariance;
    private final double lengthScale;



    public IntegratedSquaredGPApproximation(double[] coefficients,
                                            double origin,
                                            double boundary,
                                            double marginalVariance,
                                            double lengthScale) {


        this.coefficients = coefficients;
        this.origin = origin;
        this.boundary = boundary;
        this.marginalVariance = marginalVariance;
        this.lengthScale = lengthScale;

    }

    private double getConstantForIntegral(int j) {
        return Math.sqrt(getSpectralDensity(Math.sqrt(getSpectralDensityEigenValue(j, boundary)),
                marginalVariance, lengthScale)) * coefficients[j] * Math.sqrt(1/boundary);
    }


    private double getConstantForGradient(int j) {
        return Math.sqrt(getSpectralDensity(Math.sqrt(getSpectralDensityEigenValue(j, boundary)),
                marginalVariance, lengthScale)) * Math.sqrt(1/boundary);
    }


    public double getIntegral(double start, double end) {

        double length = end - start;
        double sum = origin * origin * (end -start);
        double term1;
        double term2;


        if (end > start) {

            for (int i = 0; i < coefficients.length; i++) {
                sum += getConstantForIntegral(i) * getConstantForIntegral(i) * (length/2);
                sum -= 0.25 * getConstantForIntegral(i) * getConstantForIntegral(i) * (1/Math.sqrt(getSpectralDensityEigenValue(i, boundary)))
                * (Math.sin(2 * Math.sqrt(getSpectralDensityEigenValue(i, boundary))*(end + boundary)) -
                        Math.sin(2 * Math.sqrt(getSpectralDensityEigenValue(i, boundary))*(start + boundary)));
                sum -= 2 * origin * getConstantForIntegral(i) * (1/Math.sqrt(getSpectralDensityEigenValue(i, boundary)))
                        * (Math.cos(Math.sqrt(getSpectralDensityEigenValue(i, boundary)) * (end + boundary))
                        - Math.cos(Math.sqrt(getSpectralDensityEigenValue(i, boundary)) * (start + boundary)));
            }

            for (int i = 0; i < coefficients.length; i++) {
                for (int j = i + 1; j < coefficients.length; j++) {

                    term1 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(end + boundary));
                    term2 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(start + boundary));
                    sum += (getConstantForIntegral(i) * getConstantForIntegral(j)/(Math.sqrt(getSpectralDensityEigenValue(j, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))) * (term1 - term2);
                    term1 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(end + boundary));
                    term2 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(start + boundary));
                    sum += (getConstantForIntegral(i) * getConstantForIntegral(j)/(Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))) * (term2 - term1);

                }
            }

        }

        return sum;
    }

    /*public double[] getGradient(double start, double end) {
        double[] gradients = new double[coefficients.length];
        for (int i = 0; i < coefficients.length; i++) {
            gradients[i] = getGradientWrtCoefficient(start, end, i);
        }
        return  gradients;
    }*/


    /*public double getGradientWrtCoefficient(double start, double end, int j) {

        double term1;
        double term2;
        double sum = 0;

        sum += getConstantForGradient(j) * getConstantForGradient(j) * coefficients[j]
                * ((end - start) - 0.5 * (1/Math.sqrt(getSpectralDensityEigenValue(j, boundary))) *
                (Math.sin(2 * Math.sqrt(getSpectralDensityEigenValue(j, boundary))*(end + boundary)) -
                 Math.sin(2 * Math.sqrt(getSpectralDensityEigenValue(j, boundary))*(start + boundary))));

        if (end > start) {


            for (int i = 0; i < coefficients.length; i++) {

                if (i != j) {

                    term1 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(end + boundary));
                    term2 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(start + boundary));
                    sum += ((getConstantForGradient(i) * getConstantForGradient(j))/(Math.sqrt(getSpectralDensityEigenValue(j, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))) * (term1 - term2);
                    term1 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(end + boundary));
                    term2 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(start + boundary));
                    sum += ((getConstantForGradient(i) * getConstantForGradient(j))/(Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))) * (term2 - term1);

                }
            }
        }

        return sum;

    }*/


public static void main(String[] args) {

    double[] coefficients = new double[] {1.2298 , -0.1874792 , 0.4808458 , -0.5665365 , 0.1375066,
            -1.550998 , -2.228136 , -0.232172 , 0.09171077 , -0.255557,
            0.2990763 , -0.3796941 , -0.1624707 , 0.2145575 , 0.1348038,
            0.584445 , -0.2177144 , -0.4694606 , -0.1766076 , -0.02908852,
            -0.01152857 , -0.1651964 , -0.2445626 , 0.225222 , 0.6756132,
            0.955706 , 0.1672407 , 0.004148639 , -0.1134777 , 0.4829792,
            0.03933189 , 0.667913 , 0.1893102 , 0.1800728 , 0.301086,
            -0.0201327 , 0.1331643 , 0.045907 , 0.09631207 , -0.2902463};

    double boundary = 2.5;
    double marginalVariance = 32.33321;
    double lengthScale = 0.8659804;
    double origin = 0;

    IntegratedSquaredGPApproximation approx = new IntegratedSquaredGPApproximation(coefficients, boundary, origin, marginalVariance, lengthScale);


    System.out.println("Integral 1:" + approx.getIntegral(-1, -0.5));
    System.out.println("Integral 2:" + approx.getIntegral(-0.5, 0));
    System.out.println("Integral 3:" + approx.getIntegral(0, 0.5));
    System.out.println("Integral 4:" + approx.getIntegral(0.5, 1));

}



}
