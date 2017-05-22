/*
 * ManyUniformGeoDistributionModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.geo;

import dr.geo.operators.UniformGeoSpatialOperator;
import dr.inference.distribution.AbstractDistributionLikelihood;
import dr.inference.distribution.CachedDistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.UniformOperator;
import dr.math.distributions.MultivariateDistribution;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class ManyUniformGeoDistributionModelParser extends AbstractXMLObjectParser {

    public final static String MANY_PARSER = "geoDistributionCollection";

    public String getParserName() {
        return MANY_PARSER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<GeoSpatialDistribution> distributions = new ArrayList<GeoSpatialDistribution>();
        List<Parameter> parameters = new ArrayList<Parameter>();
        List<Likelihood> oldLikelihood = new ArrayList<Likelihood>();

        for (int i = 0; i < xo.getChildCount(); ++i) {
            Object cxo = xo.getChild(i);
            MultivariateDistributionLikelihood ab = null;

            if (cxo instanceof MultivariateDistributionLikelihood) {
                ab = (MultivariateDistributionLikelihood) cxo;
            } else if (cxo instanceof CachedDistributionLikelihood) {
                ab = (MultivariateDistributionLikelihood)
                        ((CachedDistributionLikelihood) cxo).getDistributionLikelihood();
            }

            if (ab != null) {
                parameters.add(ab.getDataParameter());
                MultivariateDistribution md = ab.getDistribution();
                oldLikelihood.add((Likelihood)cxo);

                if (md instanceof GeoSpatialDistribution) {
                    distributions.add((GeoSpatialDistribution)md);
                } else {
                    throw new XMLParseException("Unhandled distribution type in '" + xo.getId() + "'");
                }

            }
        }

        return new ManyUniformGeoDistributionModel(xo.getId(), parameters, distributions, oldLikelihood);
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
            new ElementRule(MultivariateDistributionLikelihood.class, 0, Integer.MAX_VALUE),
            new ElementRule(CachedDistributionLikelihood.class, 0, Integer.MAX_VALUE),
    };
}
