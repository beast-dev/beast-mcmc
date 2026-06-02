/*
 * IndependentInverseGammaDistributionModel.java
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

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.InverseGammaDistribution;

public class IndependentInverseGammaDistributionModel extends AbstractModelLikelihood {
        Parameter shape;
        Parameter scale;
        Parameter data;

    public IndependentInverseGammaDistributionModel(String id, Parameter shape, Parameter scale, Parameter data){
        super(id);
        this.shape = shape;
        addVariable(shape);
        this.scale = scale;
        addVariable(scale);
        this.data = data;
        addVariable(data);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

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

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        double sum = 0;
        for (int i = 0; i < data.getDimension(); i++) {
            sum += InverseGammaDistribution.logPdf(data.getParameterValue(i),
                    shape.getParameterValue(i), scale.getParameterValue(i), 0.0);
        }
        return sum;
    }

    public Parameter getShape() {
        return shape;
    }

    public Parameter getScale() {
        return scale;
    }

    public Parameter getData() {
        return data;
    }

    @Override
    public void makeDirty() {

    }
}
