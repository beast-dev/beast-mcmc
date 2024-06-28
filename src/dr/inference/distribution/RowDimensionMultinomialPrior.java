/*
 * RowDimensionMultinomialPrior.java
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

/**
 * @author Max Tolkoff
 */
public class RowDimensionMultinomialPrior extends AbstractModelLikelihood implements MatrixSizePrior {
    AdaptableSizeFastMatrixParameter data;
    Parameter distribution;
    final boolean transpose;

    public RowDimensionMultinomialPrior(String id, AdaptableSizeFastMatrixParameter data, Parameter distribution, boolean transpose){
        super(id);
        addVariable(data);
        addVariable(distribution);
        this.data = data;
        this.distribution = distribution;
        this.transpose = transpose;
    }

    @Override
    public double getSizeLogLikelihood() {
        if(!transpose){
            return distribution.getParameterValue(data.getColumnDimension() - 1);
        }
        else
        {
            return distribution.getParameterValue(data.getRowDimension() - 1);
        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        return getSizeLogLikelihood();
    }

    @Override
    public void makeDirty() {
        getSizeLogLikelihood();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }
}
