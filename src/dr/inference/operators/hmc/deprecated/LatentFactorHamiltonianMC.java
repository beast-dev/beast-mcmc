/*
 * LatentFactorHamiltonianMC.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc.deprecated;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.model.*;
import dr.inference.operators.hmc.deprecated.AbstractHamiltonianMCOperator;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;

/**
 * Created by max on 12/2/15.
 */
@Deprecated
public class LatentFactorHamiltonianMC extends AbstractHamiltonianMCOperator{
    private LatentFactorModel lfm;
    private FullyConjugateMultivariateTraitLikelihood tree;
    private MatrixParameterInterface factors;
    private MatrixParameterInterface loadings;
    private MatrixParameterInterface Precision;
    private int nfac;
    private int ntaxa;
    private int ntraits;
    private double stepSize;
    private int nSteps;
    private boolean diffusionSN=true;
    private Parameter missingIndicator;


    public LatentFactorHamiltonianMC(LatentFactorModel lfm, FullyConjugateMultivariateTraitLikelihood tree, double weight, CoercionMode mode, double stepSize, int nSteps, double momentumSd){
        super(mode, momentumSd);
        setWeight(weight);
        this.lfm = lfm;
        this.tree = tree;
        this.factors = lfm.getFactors();
        this.loadings = lfm.getLoadings();
        this.Precision = lfm.getColumnPrecision();
        nfac = lfm.getFactorDimension();
        ntaxa = lfm.getFactors().getColumnDimension();
        ntraits = Precision.getRowDimension();
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.missingIndicator = lfm.getMissingIndicator();
    }



    @Override
    public double getCoercableParameter() {
        return Math.log(stepSize);
    }

    @Override
    public void setCoercableParameter(double value) {
        stepSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return stepSize;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Latent Factor Hamiltonian Monte Carlo";
    }

    private double[][] getMatrix(double[] residual){
        double answer[][] = new double[this.nfac][this.ntaxa];
        for (int i = 0; i < this.nfac ; i++) {
            for (int j = 0; j < ntraits; j++) {
                for (int k = 0; k < ntaxa; k++) {
                    if(missingIndicator == null || missingIndicator.getParameterValue(k * ntraits + j) != 1){
                        answer[i][k] -= loadings.getParameterValue(j, i) * Precision.getParameterValue(j, j) *
                            residual[k * ntraits + j];
                    }
                }
            }
        }
        return answer;
    }

    private double[][] getGradient(double[][] mean,
//                                 double[][] prec,
                                 double[] precfactor){
        double[] residual = lfm.getResidual();
        double[][] derivative = getMatrix(residual);

//        if(diffusionSN){
            for (int i = 0; i < nfac ; i++) {
                for (int j = 0; j < ntaxa; j++) {
                    derivative[i][j] += (factors.getParameterValue(i, j) - mean[j][i]) * precfactor[j];
                }

            }
//        }
//        else{
//            for (int i = 0; i <mean.length ; i++) {
//                double sumi = 0;
//                for (int j = 0; j <mean.length ; j++) {
//                    sumi += prec[i][j]*(factors.getParameterValue(j, randel) - mean[j]);
//                }
//                derivative[i] -= sumi;
//            }
//        }
        return derivative;
    }

    @Override
    public double doOperation() {



        double[][] mean = tree.getConditionalMeans();
        double precfactor[];
//        double[][] prec = null;
        double rand = MathUtils.nextDouble();
//        System.out.println(rand);
        double functionalStepSize = stepSize;

//        if(diffusionSN){
            precfactor = tree.getPrecisionFactors();
//        }
//        else {
//            prec = tree.getConditionalPrecision(randel);
//        }

        double[][] derivative = getGradient(mean,
//                prec,
                precfactor);
        drawMomentum(lfm.getFactorDimension() * ntaxa);

        double prop=0;
        for (int i = 0; i < momentum.length ; i++) {
            prop += momentum[i] * momentum[i] / (2 * getMomentumSd() * getMomentumSd());
        }


        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            for (int j = 0; j < ntaxa; j++) {
                momentum[i * ntaxa + j] = momentum[i * ntaxa + j] - functionalStepSize / 2 * derivative[i][j];
            }

        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                for (int k = 0; k < ntaxa; k++) {
                    factors.setParameterValueQuietly(j, k, factors.getParameterValue(j, k) + functionalStepSize * momentum[j * ntaxa + k] / (getMomentumSd() * getMomentumSd()));
                }
            }
            factors.fireParameterChangedEvent();

            if(i != nSteps){
                derivative = getGradient(mean, precfactor);

                for (int j = 0; j < lfm.getFactorDimension() ; j++) {
                    for (int k = 0; k < ntaxa; k++) {
                        momentum[j * ntaxa + k] = momentum[j * ntaxa + k] - functionalStepSize * derivative[j][k];
                    }
                }
            }
        }

        derivative = getGradient(mean , precfactor);
        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            for (int j = 0; j < ntaxa; j++) {
                momentum[i * ntaxa + j] = momentum[i * ntaxa + j] - functionalStepSize / 2 * derivative[i][j];
            }
        }

        double res=0;
        for (int i = 0; i <momentum.length ; i++) {
            res += momentum[i] * momentum[i] / (2 * getMomentumSd() * getMomentumSd());
        }
        return prop - res;
    }
}
