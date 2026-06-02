/*
 * ScaledParameter.java
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

package dr.inference.model;


import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Gabriel Hassler
 * @author Marc Suchard
 */


public class ScaledParameter extends ProductParameter {

    private final Parameter scaleParam;
    private final Parameter vecParam;

    public ScaledParameter(Parameter scaleParam, Parameter vecParam) {
        super(new ArrayList<>(Arrays.asList(scaleParam, vecParam)));
        this.scaleParam = scaleParam;
        this.vecParam = vecParam;

    }

    public ScaledParameter(double scale, Parameter vecParam) {
        this(new Parameter.Default(scale), vecParam);
    }

    @Override
    public double getParameterValue(int dim) {
        return scaleParam.getParameterValue(0) * vecParam.getParameterValue(dim);
    }

    @Override
    public int getDimension() {
        return vecParam.getDimension();
    }

    public Parameter getScaleParam() {
        return scaleParam;
    }

    public Parameter getVecParam() {
        return vecParam;
    }

}
