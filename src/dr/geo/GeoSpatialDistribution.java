/*
 * GeoSpatialDistribution.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateDistribution;
import dr.xml.*;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Alexei J. Drummond
 */

public class GeoSpatialDistribution implements MultivariateDistribution {

    public static final String FLAT_SPATIAL_DISTRIBUTION = "flatGeoSpatialPrior";
    public static final String DATA = "data";
    public static final String TYPE = "geoSpatial";
    public static final String NODE_LABEL = "taxon";
    public static final String KML_FILE = "kmlFileName";
    public static final String INSIDE = "inside";
    public static final String UNION = "union";
    private static final String DEFAULT_LABEL = "";

    public static final int dimPoint = 2; // Assumes 2D points only

    public GeoSpatialDistribution(String label) {
        this.label = label;
    }

    public GeoSpatialDistribution(String label, Polygon2D region, boolean inside) {
        this.label = label;
        this.region = region;
        this.outside = !inside;
    }

    public double logPdf(double[] x) {
        final boolean contains = region.containsPoint2D(new Point2D.Double(x[0], x[1]));
        if (outside ^ contains)
            return 0;
        return Double.NEGATIVE_INFINITY;

    }

    public double[][] getScaleMatrix() {
        return null;
    }

    public double[] getMean() {
        return null;
    }

    public String getType() {
        return TYPE;
    }

    public String getLabel() {
        return label;
    }

    public boolean getOutside() {
        return outside;
    }

    public Polygon2D getRegion() {
        return region;
    }

    protected Polygon2D region;
    protected String label = null;
    private boolean outside = false;

    public static XMLObjectParser FLAT_GEOSPATIAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FLAT_SPATIAL_DISTRIBUTION;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String label = xo.getAttribute(NODE_LABEL, DEFAULT_LABEL);

            boolean inside = xo.getAttribute(INSIDE, true);
            boolean union = xo.getAttribute(UNION, false);
            boolean readFromFile = false;

            List<GeoSpatialDistribution> geoSpatialDistributions = new ArrayList<GeoSpatialDistribution>();

            if (xo.hasAttribute(KML_FILE)) {
                // read file
                String kmlFileName = xo.getStringAttribute(KML_FILE);
                List<Polygon2D> polygons = Polygon2D.readKMLFile(kmlFileName);
                for (Polygon2D region : polygons)
                    geoSpatialDistributions.add(new GeoSpatialDistribution(label, region, inside));
                readFromFile = true;
            } else {

                for (int i = 0; i < xo.getChildCount(); i++) {
                    if (xo.getChild(i) instanceof Polygon2D) {
                        Polygon2D region = (Polygon2D) xo.getChild(i);
                        geoSpatialDistributions.add(
                                new GeoSpatialDistribution(label, region, inside)
                        );
                    }
                }
            }

            List<Parameter> parameters = new ArrayList<Parameter>();
            XMLObject cxo = xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                Parameter spatialParameter = (Parameter) cxo.getChild(j);
                parameters.add(spatialParameter);
            }

            if (geoSpatialDistributions.size() == 1) {
                MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(geoSpatialDistributions.get(0));
                for (Parameter spatialParameter : parameters) {
                    if (spatialParameter.getDimension() != dimPoint)
                        throw new XMLParseException("Spatial priors currently only work in " + dimPoint + "D");
                    likelihood.addData(spatialParameter);
                }
                return likelihood;
            }

            if (geoSpatialDistributions.size() == 0) {
                throw new XMLParseException("Error constructing geo spatial distributions in " + xo.getId());
            }

            if (parameters.size() == 1) {
                Parameter parameter = parameters.get(0);
                if (parameter.getDimension() % dimPoint != 0)
                    throw new XMLParseException("Spatial priors currently only work in " + dimPoint + "D");

                if (!label.equals(DEFAULT_LABEL)) {  // For a tip-taxon
                    Logger.getLogger("dr.geo").info(
                            "\nConstructing a multiple-region spatial prior:\n" +
                                    "\tTaxon: " + label + "\n" +
                                    "\tNumber of regions: " + geoSpatialDistributions.size() + "\n\n");
                    MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
                            new MultiRegionGeoSpatialDistribution(label, geoSpatialDistributions, union));
                    likelihood.addData(parameter);
                    return likelihood;

                } else {

                    Logger.getLogger("dr.geo").info(
                            "\nConstructing a GeoSpatialCollectionModel:\n" +
                                    "\tParameter: " + parameter.getId() + "\n" +
                                    "\tNumber of regions: " + geoSpatialDistributions.size() + "\n\n");

                    return new GeoSpatialCollectionModel(xo.getId(), parameter, geoSpatialDistributions, !union);
                }
            }

            throw new XMLParseException("Multiple separate parameters and multiple regions not yet implemented");
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(NODE_LABEL, true),
                AttributeRule.newBooleanRule(INSIDE, true),
                AttributeRule.newBooleanRule(UNION, true),
                new XORRule(
                        AttributeRule.newStringRule(KML_FILE),
                        new ElementRule(Polygon2D.class, 1, Integer.MAX_VALUE)
                ),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)}
                )
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under a 2D geospatial distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

}
