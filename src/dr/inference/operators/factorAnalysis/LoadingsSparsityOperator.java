/*
 * LoadingsSparsityOperator.java
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

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.factorAnalysis.LoadingsGibbsTruncatedOperator;
import dr.math.MathUtils;

public class LoadingsSparsityOperator extends SimpleMCMCOperator {
    public LoadingsSparsityOperator(double weight, LoadingsGibbsTruncatedOperator loadingsGibbs, MatrixParameterInterface sparse){

        setWeight(weight);
        this.loadingsGibbs = loadingsGibbs;
        this.sparse = sparse;
    }

    @Override
    public String getOperatorName() {
        return "LoadingsSparsityOperator";
    }

    @Override
    public double doOperation() {
        int row = MathUtils.nextInt(sparse.getRowDimension());
        int col = MathUtils.nextInt(sparse.getColumnDimension());

        double hastings = 0;
        hastings += loadingsGibbs.drawI(row, col, false);
        if(sparse.getParameterValue(row, col) == 0)
            sparse.setParameterValue(row, col, 1);
        else
            sparse.setParameterValue(row, col, 0);
        hastings -= loadingsGibbs.drawI(row, col, true);

        if(Double.isNaN(hastings)){
            System.out.println("is NaN");
        }
        if(Double.isInfinite(hastings))
        {
            return Double.NEGATIVE_INFINITY;
        }
        return hastings;
    }


    LoadingsGibbsTruncatedOperator loadingsGibbs;
    MatrixParameterInterface sparse;
}
