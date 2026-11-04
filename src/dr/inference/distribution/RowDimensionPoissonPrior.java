/*
 * RowDimensionPoissonPrior.java
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

package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.distributions.PoissonDistribution;

/**
 * Created by maxryandolinskytolkoff on 7/20/16.
 */
public class RowDimensionPoissonPrior extends AbstractModelLikelihood implements MatrixSizePrior {
    public RowDimensionPoissonPrior(String id, double untruncatedMean, MatrixParameterInterface parameter, DeterminentalPointProcessPrior DPP, boolean transpose){
        super(id);
        this.poisson = new PoissonDistribution(untruncatedMean);
        this.parameter = parameter;
        if(parameter != null)
            addVariable(parameter);
        this.transpose = transpose;
        this.DPP = DPP;
        if(parameter != null && !(parameter instanceof AdaptableSizeFastMatrixParameter)){
            if(!transpose){
                count = new int[parameter.getColumnDimension()];
                storedCount = new int[parameter.getColumnDimension()];
                for (int i = 0; i < count.length; i++) {
                    for (int j = 0; j < parameter.getRowDimension(); j++) {
                        count[i] += parameter.getParameterValue(j, i);
                    }

                }
            }
            else {
                count = new int[parameter.getRowDimension()];
                storedCount = new int[parameter.getRowDimension()];
                for (int i = 0; i < count.length; i++) {
                    for (int j = 0; j < parameter.getColumnDimension(); j++) {
                        count[i] += parameter.getParameterValue(i, j);
                    }

                }
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState(){
        if(count != null)
            System.arraycopy(count, 0, storedCount, 0, count.length);
    }

    @Override
    protected void restoreState() {
        int[] temp = storedCount;
        storedCount = count;
        count = temp;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        System.out.println("Variable Changed");
//        System.out.println(index);
//        if(parameter != null && !(parameter instanceof AdaptableSizeFastMatrixParameter)){
//            int row;
//            int col;
//            if(!transpose){
//                row = index % parameter.getRowDimension();
//                col = index / parameter.getRowDimension();
//                if(parameter.getParameterValue(row, col) == 1){
//                    count[col] += 1;
//                }
//                else{
//                    count[col] -= 1;
//                    System.out.println(row);
//                    System.out.println(col);
//                    System.out.println(count[col]);
//                    System.out.println(parameter.getParameterValue(row, col));
//                }
//            }
//            else{
//                if(!transpose){
//                    row = index / parameter.getColumnDimension();
//                    col = index % parameter.getColumnDimension();
//                    if(parameter.getParameterValue(row, col) == 1){
//                        count[row] += 1;
//                    }
//                    else
//                        count[row] -= 1;
//                }
//            }
//        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
//        System.out.println(poisson.logPdf(parameter.getRowDimension()) - Math.log(1 - Math.exp(-poisson.mean())));
        if (DPP != null)
            return poisson.logPdf(DPP.getSum()) - Math.log(1 - Math.exp(-poisson.mean()));
        if(count != null){
            if(parameter != null && !(parameter instanceof AdaptableSizeFastMatrixParameter)){
                if(!transpose){
                    count = new int[parameter.getColumnDimension()];
                    storedCount = new int[parameter.getColumnDimension()];
                    for (int i = 0; i < count.length; i++) {
                        for (int j = 0; j < parameter.getRowDimension(); j++) {
                            count[i] += parameter.getParameterValue(j, i);
                        }

                    }
                }
                else {
                    count = new int[parameter.getRowDimension()];
                    storedCount = new int[parameter.getRowDimension()];
                    for (int i = 0; i < count.length; i++) {
                        for (int j = 0; j < parameter.getColumnDimension(); j++) {
                            count[i] += parameter.getParameterValue(i, j);
                        }

                    }
                }
            }

            double nonZero = 0;
            for (int i = 0; i < count.length; i++) {
                if(count[i] != 0){
                    nonZero ++;
                }
                if(count[i] < 0)
                    throw new RuntimeException("Less than 0 elements are 0. Check for an error");
            }
            return poisson.logPdf(nonZero) - Math.log(1 - Math.exp(-poisson.mean()));
        }
        if(!transpose){
            return poisson.logPdf(parameter.getRowDimension()) - Math.log(1 - Math.exp(-poisson.mean()));
        }
        else
            return poisson.logPdf(parameter.getColumnDimension()) - Math.log(1 - Math.exp(-poisson.mean()));
    }

    @Override
    public void makeDirty() {

    }

    @Override
    public double getSizeLogLikelihood() {
        return getLogLikelihood();
    }

    PoissonDistribution poisson;
    MatrixParameterInterface parameter;
    boolean transpose;
    DeterminentalPointProcessPrior DPP;
    int[] count;
    int[] storedCount;
}
