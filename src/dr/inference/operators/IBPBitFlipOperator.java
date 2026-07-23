/*
 * IBPBitFlipOperator.java
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

package dr.inference.operators;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.IndianBuffetProcessPrior;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * Created by maxryandolinskytolkoff on 11/21/16.
 */
public class IBPBitFlipOperator extends BitFlipOperator {
    public IBPBitFlipOperator(Parameter parameter, double weight, IndianBuffetProcessPrior IBP) {
        super(parameter, weight, true);
        this.sparsity = (MatrixParameterInterface) parameter;
        this.IBP = IBP;
    }

    @Override
    public double doOperation(Likelihood likelihood)
    { //throws OperatorFailedException {
        return super.doOperation(likelihood);
    }

    protected double sum(int pos) {
        double sum = 0;
        int column = pos / sparsity.getRowDimension();
        for (int i = 0; i < sparsity.getRowDimension(); i++) {
            sum += sparsity.getParameterValue(i, column);
        }


        return sum;
    }

    protected  int getDimension(){
        return sparsity.getRowDimension();
    }

    IndianBuffetProcessPrior IBP;
    MatrixParameterInterface sparsity;


}
