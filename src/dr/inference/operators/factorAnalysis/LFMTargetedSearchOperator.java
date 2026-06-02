/*
 * LFMTargetedSearchOperator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.operators.factorAnalysis;

import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;


import java.util.ArrayList;

public class LFMTargetedSearchOperator extends SimpleMCMCOperator {
    public LFMTargetedSearchOperator(Double weight, MatrixParameterInterface SparseMatrix, ArrayList<MatrixParameterInterface> MatrixList,
    MatrixParameterInterface FactorsMatrix, ArrayList<MatrixParameterInterface> FactorsMatrixList,
    MatrixParameterInterface LoadingsMatrix, ArrayList<MatrixParameterInterface> LoadingsMatrixList,
                                     MatrixParameterInterface cutoffs){
        setWeight(weight);
        this.SparseMatrix = SparseMatrix;
        this.MatrixList = MatrixList;
        this.FactorsMatrix = FactorsMatrix;
        this.FactorsMatrixList = FactorsMatrixList;
        this.LoadingsMatrix = LoadingsMatrix;
        this.LoadingsMatrixList = LoadingsMatrixList;
//        Parameter mean = new Parameter.Default(1);
//        mean.setParameterValue(0, 0);
//        Parameter SD = new Parameter.Default(1);
//        SD.setParameterValue(0, this.VAR);
        this.error = new NormalDistribution(0, VAR);
        this.cutoffs = cutoffs;
    }

    @Override
    public String getOperatorName() {
        return "Latent Factor Model Targeted Search Operator";
    }

    @Override
    public double doOperation() {
        final int DRAW = MathUtils.nextInt(MatrixList.size());
        lastDraw = DRAW;
        double hastings = getHastings();
        for (int i = 0; i < SparseMatrix.getDimension(); i++) {
            if(MathUtils.nextDouble() < PROBABILITY)
                SparseMatrix.setParameterValue(i, MatrixList.get(DRAW).getParameterValue(i));
            else
                SparseMatrix.setParameterValue(i, MathUtils.nextInt(2));
        }
        for (int i = 0; i < FactorsMatrix.getDimension(); i++) {
            FactorsMatrix.setParameterValue(i, FactorsMatrixList.get(DRAW).getParameterValue(i) + (Double) error.nextRandom());
        }
        for (int i = 0; i < LoadingsMatrix.getDimension(); i++) {
            boolean badDraw = true;
            int count = 0;
            while(badDraw || count == 10){
                double draw = LoadingsMatrixList.get(DRAW).getParameterValue(i) + (Double) error.nextRandom();
                if(-Math.sqrt(cutoffs.getParameterValue(i)) < draw && Math.sqrt(cutoffs.getParameterValue(i)) > draw) {
                    badDraw = false;
                    LoadingsMatrix.setParameterValue(i, draw);
                }
                else{
                    count ++;
                }
            }
        }

        hastings -= getHastings();


        return hastings;
    }

    private double getHastings(){
        int denom = MatrixList.size();
        double hastings = 0;
        for (int i = 0; i < SparseMatrix.getDimension(); i++) {
            int sum = 0;
            for (int j = 0; j < MatrixList.size(); j++) {
                if(MatrixList.get(j).getParameterValue(i) != SparseMatrix.getParameterValue(i)){
                    sum ++;
                }
            }
            hastings += Math.log((sum / denom) * PROBABILITY + .5 * (1 - PROBABILITY));
        }
        for (int i = 0; i < FactorsMatrix.getDimension(); i++) {
            double sum = 0;
            for (int j = 0; j < FactorsMatrixList.size(); j++) {
                sum += error.logPdf(FactorsMatrix.getParameterValue(i) - FactorsMatrixList.get(j).getParameterValue(i)) - Math.log(FactorsMatrixList.size());
            }
            hastings += sum;
        }
        for (int i = 0; i < LoadingsMatrix.getDimension(); i++) {
            double sum = 0;
            for (int j = 0; j < LoadingsMatrixList.size(); j++) {
                sum += error.logPdf(LoadingsMatrix.getParameterValue(i) - LoadingsMatrixList.get(j).getParameterValue(i)) -
                        Math.log(error.cdf(Math.sqrt(cutoffs.getParameterValue(i)) - LoadingsMatrixList.get(j).getParameterValue(i))
                                - error.cdf(-Math.sqrt(cutoffs.getParameterValue(i)) - LoadingsMatrixList.get(j).getParameterValue(i)))
                        - Math.log(FactorsMatrixList.size());
            }
            hastings += sum;
        }


        return hastings;
    }


    @Override
    public void accept(double deviation) {
        super.accept(deviation);
        System.out.println(lastDraw);
    }

    ArrayList<MatrixParameterInterface> MatrixList;
    MatrixParameterInterface SparseMatrix;
    ArrayList<MatrixParameterInterface> FactorsMatrixList;
    MatrixParameterInterface FactorsMatrix;
    ArrayList<MatrixParameterInterface> LoadingsMatrixList;
    MatrixParameterInterface LoadingsMatrix;
    NormalDistribution error;
    MatrixParameterInterface cutoffs;
    private final double PROBABILITY = .9;
    private final double VAR = .15;
    private int lastDraw;
}
