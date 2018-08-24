/*
 * LatentFactorModelPrecisionGibbsOperator.java
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

package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;

/**
 * Created by max on 6/12/14.
 */
public class LatentFactorModelPrecisionGibbsOperator extends SimpleMCMCOperator implements PathDependentOperator, GibbsOperator {
    //    private double[] FacXLoad;
//    private double[] residual;
    private LatentFactorModel LFM;
    private GammaDistribution prior;
    private boolean randomScan;
    private double shape;
    double pathWeight = 1.0;
    final Parameter missingIndicator;

    public LatentFactorModelPrecisionGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight, boolean randomScan) {
        setWeight(weight);
        this.LFM = LFM;
        this.prior = (GammaDistribution) prior.getDistribution();
        this.randomScan = randomScan;
        this.missingIndicator = LFM.getMissingIndicator();

//        FacXLoad=new double[LFM.getFactors().getColumnDimension()];
//        residual=new double[LFM.getFactors().getColumnDimension()];
//        setShape();
    }

    private double getShape(int i) {
        return this.prior.getShape() + LFM.getRowCount(i) * .5 * pathWeight;
    }

    private void setPrecision(int i) {
        MatrixParameterInterface factors = LFM.getFactors();
        MatrixParameterInterface loadings = LFM.getLoadings();
        DiagonalMatrix precision = (DiagonalMatrix) LFM.getColumnPrecision();
        MatrixParameterInterface data = LFM.getScaledData();
        double di = 0;
        for (int j = 0; j < factors.getColumnDimension(); j++) {
            double sum = 0;
            if(missingIndicator == null || missingIndicator.getParameterValue(j * LFM.getScaledData().getRowDimension() + i) != 1) {
                for (int k = 0; k < factors.getRowDimension(); k++) {
                    sum += factors.getParameterValue(k, j) * loadings.getParameterValue(i, k);
                }
                double temp = data.getParameterValue(i, j) - sum;
                di += temp * temp;
            }
//            FacXLoad[j]=sum;
//            residual[j]=data.getParameterValue(i,j)-FacXLoad[j];
//        }
//        double sum=0;
//        for (int j = 0; j <factors.getColumnDimension() ; j++) {
//            sum+=residual[j]*residual[j];
        }
        double scale = 1.0 / (1.0 / prior.getScale() + pathWeight * di * .5);
        double nextPrecision = GammaDistribution.nextGamma(getShape(i), scale);
        precision.setParameterValueQuietly(i, nextPrecision);
    }

    public void setPathParameter(double beta) {
        pathWeight = beta;
    }

    @Override
    public String getPerformanceSuggestion() {
        return "Only works for diagonal column precision matrices for a LatentFactorModel with a gamma prior";
    }

    @Override
    public String getOperatorName() {
        return "Latent Factor Model Precision Gibbs Operator";
    }

    @Override
    public double doOperation() {

        if (!randomScan) for (int i = 0; i < LFM.getColumnPrecision().getColumnDimension(); i++) {
            if (LFM.getContinuous().getParameterValue(i) != 0)
                setPrecision(i);
        }
        else {
            int i = MathUtils.nextInt(LFM.getColumnPrecision().getColumnDimension());
            while (LFM.getContinuous().getParameterValue(i) == 0)
                i = MathUtils.nextInt(LFM.getColumnPrecision().getColumnDimension());
            setPrecision(i);
        }
        LFM.getColumnPrecision().getParameter(0).fireParameterChangedEvent();


        return 0;
    }
}
