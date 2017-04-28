/*
 * LoadingsHamiltonianMC.java
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

import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.*;
import dr.inference.operators.hmc.deprecated.AbstractHamiltonianMCOperator;
import dr.inference.operators.CoercionMode;

/**
 * Created by max on 1/11/16.
 */
@Deprecated
public class LoadingsHamiltonianMC extends AbstractHamiltonianMCOperator {
    private LatentFactorModel lfm;
    private MomentDistributionModel prior;
    private MatrixParameterInterface factors;
    private MatrixParameterInterface loadings;
    private MatrixParameterInterface Precision;
    private int nfac;
    private int ntaxa;
    private int ntraits;
    private double stepSize;
    private int nSteps;


    public LoadingsHamiltonianMC(LatentFactorModel lfm, MomentDistributionModel prior, double weight, CoercionMode mode, double stepSize, int nSteps, double momentumSd, MatrixParameterInterface loadings){
        super(mode , momentumSd);
        setWeight(weight);
        this.lfm = lfm;
        this.prior = prior;
        this.factors = lfm.getFactors();
        this.loadings = loadings;
        this.Precision = lfm.getColumnPrecision();
        nfac = lfm.getFactorDimension();
        ntaxa = lfm.getFactors().getColumnDimension();
        ntraits = Precision.getRowDimension();
        this.stepSize = stepSize;
        this.nSteps = nSteps;
    }

    @Override
    public double getCoercableParameter() {
        return 0;
    }

    @Override
    public void setCoercableParameter(double value) {

    }

    @Override
    public double getRawParameter() {
        return 0;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "LoadingsHamiltonianMC";
    }

    @Override
    public double doOperation() {


        double[][] derivative = getGradient();
        drawMomentum(lfm.getFactorDimension() * ntraits);
        double functionalStepSize = stepSize;

        double prop=0;
        for (int i = 0; i < momentum.length ; i++) {
            prop += momentum[i] * momentum[i] / (2 * getMomentumSd() * getMomentumSd());
        }

        for (int i = 0; i < lfm.getFactorDimension() ; i++) {
            for (int j = 0; j < ntraits ; j++) {
                momentum[i * ntraits + j] = momentum[i * ntraits + j] - functionalStepSize / 2 * derivative[j][i];
            }

        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                for (int k = 0; k <ntraits ; k++) {
                    loadings.setParameterValueQuietly(k, j, loadings.getParameterValue(k, j) + functionalStepSize * momentum[j * ntraits + k]);                }

            }
            loadings.fireParameterChangedEvent(-1, Parameter.ChangeType.ALL_VALUES_CHANGED);


            if(i != nSteps){
                derivative = getGradient();

                for (int j = 0; j < lfm.getFactorDimension() ; j++) {
                    for (int k = 0; k <ntraits ; k++) {
                        momentum[j * ntraits + k] = momentum[j * ntraits + k] - functionalStepSize * derivative[k][j];
                    }

                }
            }
        }

        derivative=getGradient();
        for (int i = 0; i < lfm.getFactorDimension() ; i++) {
            for (int j = 0; j < ntraits ; j++) {
                momentum[i * ntraits + j] = momentum[i * ntraits + j ] - functionalStepSize / 2 * derivative[j][i];
            }

        }
        double res=0;
        for (int i = 0; i <momentum.length ; i++) {
            res+=momentum[i] * momentum[i] / (2 * getMomentumSd() * getMomentumSd());
        }
        return prop - res;
    }



    private double[][] getLFMDerivative(){
        double[] residual=lfm.getResidual();
        double[][] answer= new double[ntraits][lfm.getFactorDimension()];
        for (int i = 0; i < ntaxa; i++) {
            for (int j = 0; j < ntraits; j++) {
                for (int k = 0; k < lfm.getFactorDimension() ; k++) {
                    answer[j][k] -= residual[i * ntaxa + j] * factors.getParameterValue(k , i);
                }
            }

        }
        for (int i = 0; i < ntraits ; i++) {
            for (int j = 0; j < lfm.getFactorDimension() ; j++) {
                answer[i][j] *= Precision.getParameterValue(i , i);
            }

        }
        return answer;
    }

    private double[][] getGradient(){
        double[][] answer = getLFMDerivative();
        for (int i = 0; i < loadings.getRowDimension(); i++) {
            for (int j = 0; j < loadings.getColumnDimension(); j++) {
                answer[i][j] += 2 / loadings.getParameterValue(i, j) + (loadings.getParameterValue(i, j) - prior.getMean()[0]) / prior.getScaleMatrix()[0][0];
            }

        }
        return answer;
    }

}
