/*
 * UniformGeoSpatialOperatorParser.java
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

package dr.geo.operators;

import dr.geo.AbstractPolygon2D;
import dr.geo.GeoSpatialDistribution;
import dr.geo.MultiRegionGeoSpatialDistribution;
import dr.geo.Polygon2D;
import dr.inference.distribution.AbstractDistributionLikelihood;
import dr.inference.distribution.CachedDistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.UniformOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class UniformGeoSpatialOperatorParser extends AbstractXMLObjectParser {
    public final static String UNIFORM_OPERATOR = "uniformGeoSpatialOperator";
    public static final String LOWER = "lower";
    public static final String UPPER = "upper";

    public String getParserName() {
        return UNIFORM_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        if( parameter.getDimension() == 0 ) {
             throw new XMLParseException("parameter with 0 dimension.");
        }



        MultivariateDistributionLikelihood likelihood = (MultivariateDistributionLikelihood)
                xo.getChild(MultivariateDistributionLikelihood.class);

        if (likelihood == null) {
            CachedDistributionLikelihood cached = (CachedDistributionLikelihood) xo.getChild(CachedDistributionLikelihood.class);
            AbstractDistributionLikelihood ab = cached.getDistributionLikelihood();
            if (!(ab instanceof MultivariateDistributionLikelihood)) {
                throw new XMLParseException("invalid likelihood type in " + xo.getId());
            }

            likelihood = (MultivariateDistributionLikelihood) ab;
        }

        List<AbstractPolygon2D> polygonList = new ArrayList<AbstractPolygon2D>();

        if (likelihood.getDistribution() instanceof MultiRegionGeoSpatialDistribution) {
            for (GeoSpatialDistribution spatial : ((MultiRegionGeoSpatialDistribution) likelihood.getDistribution()).getRegions()) {
                polygonList.add(spatial.getRegion());
            }
        } else if (likelihood.getDistribution() instanceof GeoSpatialDistribution) {
            polygonList.add(
                    ((GeoSpatialDistribution) likelihood.getDistribution()).getRegion()
            );
        } else {
            throw new XMLParseException("Multivariate distribution must be either a GeoSpatialDistribution " +
                "or a MultiRegionGeoSpatialDistribution");
        }

        return new UniformGeoSpatialOperator(parameter, weight, polygonList);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "An operator that picks new parameter values uniformly at random.";
    }

    public Class getReturnType() {
        return UniformOperator.class;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
//            AttributeRule.newDoubleRule(LOWER, true),
//            AttributeRule.newDoubleRule(UPPER, true),
            new ElementRule(Parameter.class),
            new XORRule(
            new ElementRule(MultivariateDistributionLikelihood.class),
                    new ElementRule(CachedDistributionLikelihood.class)
                    ),
    };
}
