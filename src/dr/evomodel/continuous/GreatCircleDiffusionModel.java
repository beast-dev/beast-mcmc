/*
 * GreatCircleDiffusionModel.java
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

package dr.evomodel.continuous;

import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class GreatCircleDiffusionModel extends MultivariateDiffusionModel {

    public static final String DIFFUSION_PROCESS = "greatCircleDiffusionModel";
    public static final String DIFFUSION_CONSTANT = "precision";
    public static final String COEFFICIENT = "diffusionCoefficient";
//    public static final String BIAS = "mu";
//    public static final String PRECISION_TREE_ATTRIBUTE = "precision";

    public GreatCircleDiffusionModel(Parameter precision, Parameter coefficient) {
        super();
        this.precision = precision;
        addVariable(precision);
        this.coefficient = coefficient;
        if (coefficient != null)
            addVariable(coefficient);
    }

    public GreatCircleDiffusionModel(Parameter precision) {
        this(precision, null);
    }

    protected double calculateLogDensity(double[] start, double[] stop, double time) {

        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(start[0], start[1]);
        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(stop[0], stop[1]);

        double distance = coord1.distance(coord2);

        double inverseVariance = precision.getParameterValue(0) / time;
        // TODO Check!
        // I believe this is a 2D (not 1D) Normal diffusion approx; hence the precision is squared (compared to 1D) 
        // in the normalization constant
        if (coefficient == null)
            return -LOG2PI + Math.log(inverseVariance) - 0.5 * (distance * distance * inverseVariance);
        double coef = -coefficient.getParameterValue(0);
        return -LOG2PI + coef * Math.log(inverseVariance) - 0.5 * (distance * distance * Math.pow(inverseVariance, coef));
    }

    protected void calculatePrecisionInfo() {
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter diffusionParam = (Parameter) xo.getChild(Parameter.class);
            Parameter coefficient = null;
            if (xo.hasChildNamed(COEFFICIENT))
                coefficient = (Parameter) xo.getChild(COEFFICIENT).getChild(Parameter.class);

            return new GreatCircleDiffusionModel(diffusionParam, coefficient);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Describes a bivariate diffusion process using great circle distances.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class),
                new ElementRule(COEFFICIENT,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
        };

        public Class getReturnType() {
            return MultivariateDiffusionModel.class;
        }
    };

    private Parameter precision;
    private Parameter coefficient;

}
