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
    private final double degree;
    private final double boundary;
    private final double mean;
    private final double marginalVariance;
    private final double lengthScale;


    private double integral;
    private double storedIntegral;
    private boolean storedIntegralKnown;


    public IntegratedSquaredGPApproximation(double[] coefficients,
                                            double degree,
                                            double boundary,
                                            double mean,
                                            double marginalVariance,
                                            double lengthScale) {


        this.coefficients = coefficients;
        this.boundary = boundary;
        this.mean = mean;
        this.marginalVariance = marginalVariance;
        this.lengthScale = lengthScale;
        this.degree = degree;
    }

    private double getConstantForIntegral(int j) {
        return Math.sqrt(getSpectralDensity(Math.sqrt(getSpectralDensityEigenValue(j, boundary)),
                marginalVariance, lengthScale, degree)) * coefficients[j] * Math.sqrt(1/boundary);
    }


    private double getConstantForGradient(int j) {
        return Math.sqrt(getSpectralDensity(Math.sqrt(getSpectralDensityEigenValue(j, boundary)),
                marginalVariance, lengthScale, degree)) * Math.sqrt(1/boundary);
    }


    public double getIntegral(double start, double end) {

        double length = end - start;
        double sum = mean * mean * (end - start);
        double term1;
        double term2;


        if (end > start) {

            for (int i = 0; i < coefficients.length; i++) {
                sum += getConstantForIntegral(i) * getConstantForIntegral(i) * (length/2);
                sum -= 0.25 * getConstantForIntegral(i) * getConstantForIntegral(i) *
                        (1/Math.sqrt(getSpectralDensityEigenValue(i, boundary)))
                        * (Math.sin(2 * Math.sqrt(getSpectralDensityEigenValue(i, boundary))*(end + boundary)) -
                        Math.sin(2 * Math.sqrt(getSpectralDensityEigenValue(i, boundary))*(start + boundary)));
                sum -= 2 * mean * getConstantForIntegral(i) * (1/Math.sqrt(getSpectralDensityEigenValue(i, boundary))) *
                        (Math.cos(Math.sqrt(getSpectralDensityEigenValue(i, boundary)) * (end + boundary))
                        -Math.cos(Math.sqrt(getSpectralDensityEigenValue(i, boundary)) * (start + boundary)));
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



    public double getGradientWrtCoefficient(double start, double end, int j) {

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

                    term1 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(i, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(j, boundary)))*(end + boundary));
                    term2 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(i, boundary)) -
                            Math.sqrt(getSpectralDensityEigenValue(j, boundary)))*(start + boundary));
                    sum += ((getConstantForGradient(i) * getConstantForGradient(j) * coefficients[i])/
                            (Math.sqrt(getSpectralDensityEigenValue(i, boundary)) -
                                    Math.sqrt(getSpectralDensityEigenValue(j, boundary)))) * (term1 - term2);
                    term1 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(end + boundary));
                    term2 = Math.sin((Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                            Math.sqrt(getSpectralDensityEigenValue(i, boundary)))*(start + boundary));
                    sum += ((getConstantForGradient(i) * getConstantForGradient(j) * coefficients[i])/
                            (Math.sqrt(getSpectralDensityEigenValue(j, boundary)) +
                                    Math.sqrt(getSpectralDensityEigenValue(i, boundary)))) * (term2 - term1);

                }
            }
        }

        return sum;

    }

    public static void main(String[] args) {

        double[] coefficients = new double[] {1.686465, -47.05261, 2.638783, -9.978589, 1.249006,
                8.846663, 2.036995, 9.899196, 0.3506009, 3.975602,
                -1.652094, 2.907649, -7.169512, -1.804519, -4.756368,
                -28.3502, 23.20718, 28.51416, -5.019321, 4.746734};

        double boundary = 1.5;
        double mean = 19.52;
        double marginalVariance = 1.0;
        double lengthScale = 0.039;
        double degree = 0.5;

        IntegratedSquaredGPApproximation approx = new IntegratedSquaredGPApproximation(coefficients, degree, boundary,
                mean, marginalVariance, lengthScale);


        System.out.println("Integral 1:" + approx.getIntegral(-1, -0.5));
        System.out.println("Integral 2:" + approx.getIntegral(-0.5, 0));
        System.out.println("Integral 3:" + approx.getIntegral(0, 0.5));
        System.out.println("Integral 4:" + approx.getIntegral(0.5, 1));

    }



}