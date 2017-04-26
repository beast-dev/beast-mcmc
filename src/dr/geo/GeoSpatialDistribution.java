/*
 * GeoSpatialDistribution.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.distribution.CachedDistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
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

public class GeoSpatialDistribution implements MultivariateDistribution, Citable {

    public static final String FLAT_SPATIAL_DISTRIBUTION = "flatGeoSpatialPrior";
    public static final String DATA = "data";
    public static final String TYPE = "geoSpatial";
    public static final String NODE_LABEL = "taxon";
    public static final String KML_FILE = "kmlFileName";
    public static final String INSIDE = "inside";
    public static final String UNION = "union";
    public static final String CACHE = "cache";
    private static final String DEFAULT_LABEL = "";

    public static final int dimPoint = 2; // Assumes 2D points only

    public GeoSpatialDistribution(String label) {
        this.label = label;
    }

    public GeoSpatialDistribution(String label, AbstractPolygon2D region, boolean inside) {
        this.label = label;
        this.region = region;
        this.outside = !inside;
    }

    public double logPdf(double[] x) {
        /*final boolean contains = region.containsPoint2D(new Point2D.Double(x[0], x[1]));
        if (region.hasFillValue()) {
            if (contains) {
                return region.getLogFillValue();
            } else {
                return Double.NEGATIVE_INFINITY;
            }
        } else if (outside ^ contains)
            return 0;
        return Double.NEGATIVE_INFINITY;*/

        return region.getLogProbability(new Point2D.Double(x[0], x[1]), this.outside);
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

    public AbstractPolygon2D getRegion() {
        return region;
    }

    protected AbstractPolygon2D region;
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
            boolean cache = xo.getAttribute(CACHE, false);

            List<GeoSpatialDistribution> geoSpatialDistributions = new ArrayList<GeoSpatialDistribution>();

            boolean overAllFill = false;
            if (xo.hasAttribute(KML_FILE)) {
                // read file
                String kmlFileName = xo.getStringAttribute(KML_FILE);
                List<AbstractPolygon2D> polygons = Polygon2D.readKMLFile(kmlFileName);
                //check if all the polygons have either a fillValue or do not have this fillValue
                AbstractPolygon2D first = polygons.get(0);
                boolean equalFills = first.hasFillValue();
                if (polygons.size() > 1) {
                    for (int i = 1; i < polygons.size(); i++) {
                        if (!(equalFills == polygons.get(i).hasFillValue())) {
                            throw new XMLParseException("Inconsistent fillValue attributes provided.");
                        }
                    }
                }
                if (equalFills) {
                    double sum = 0.0;
                    for (AbstractPolygon2D region : polygons) {
                        sum += region.getFillValue();
                    }
                    if (Math.abs(sum - 1.0) > 1E-12) {
                        throw new XMLParseException("Fill values in " + kmlFileName + " do not sum to 1 : " + sum);
                    }
                }
                overAllFill = equalFills;
                for (AbstractPolygon2D region : polygons) {
                    geoSpatialDistributions.add(new GeoSpatialDistribution(label, region, inside));
                }
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

            //changes with fillValue attribute do not apply here
            if (geoSpatialDistributions.size() == 1) {
                MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(geoSpatialDistributions.get(0));
                for (Parameter spatialParameter : parameters) {
                    if (spatialParameter.getDimension() % dimPoint != 0)
                        throw new XMLParseException("Spatial priors currently only work in " + dimPoint + "D");

                    if (!label.equals(DEFAULT_LABEL)) {  // For a tip-taxon
                        likelihood.addData(spatialParameter);
                        return likelihood;
                    } else {
                        Logger.getLogger("dr.geo").info(
                                "\nConstructing a GeoSpatialCollectionModel:\n" +
                                        "\tParameter: " + spatialParameter.getId() + "\n" +
                                        "\tNumber of regions: " + geoSpatialDistributions.size() + "\n\n");

                        return new GeoSpatialCollectionModel(xo.getId(), spatialParameter, geoSpatialDistributions, !union);
                    }
                }
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
                            new MultiRegionGeoSpatialDistribution(label, geoSpatialDistributions, union, overAllFill));
                    likelihood.addData(parameter);
                    likelihood.setId(xo.getId());
                    if (cache) {
                        return new CachedDistributionLikelihood(xo.getId(), likelihood, parameter);
                    } else {
                        return likelihood;
                    }

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
                AttributeRule.newBooleanRule(CACHE, true),
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
            return MultivariateDistributionLikelihood.class;
        }
    };

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "Integrated continuous traits over polygons";
    }

    public List<Citation> getCitations() {
        List<Citation> citationList = new ArrayList<Citation>();
        citationList.add(new Citation(
                new Author[] {
                        new Author("S", "Nylinder"),
                        new Author("P", "Lemey"),
                        new Author("M", "de Bruyn"),
                        new Author("MA", "Suchard"),
                        new Author("BE", "Pfeil"),
                        new Author("N", "Walsh"),
                        new Author("AA", "Anderberg")
                },
                "On the biogeography of Centipeda: a species-tree diffusion approach",
                2014,
                "Systematic Biology",
                63,
                178, 191,
                Citation.Status.PUBLISHED
        ));
        return citationList;
    }

    public static void main(String[] args) {

        //test reading in KML file with multiple polygons and fillValue attribute
        List<AbstractPolygon2D> polygons = Polygon2D.readKMLFile("Multiple_polygons/9.285_105.7244444.kml");
        System.out.println(polygons.size() + " polygons found.");

        //check if all the polygons have either a fillValue or do not have this fillValue
        AbstractPolygon2D first = polygons.get(0);
        boolean equalFills = first.hasFillValue();
        if (polygons.size() > 1) {
            for (int i = 1; i < polygons.size(); i++) {
                if (!(equalFills == polygons.get(i).hasFillValue())) {
                    System.out.println("Inconsistent fillValue attributes provided.");
                }
            }
        }

        if (equalFills) {
            double sum = 0.0;
            for (AbstractPolygon2D region : polygons) {
                sum += region.getFillValue();
            }
            if (Math.abs(sum - 1.0) > 1E-12) {
                System.out.println("\nFill values do not sum to 1 : " + sum);
            }
        }

        boolean overAllFill = equalFills;
        List<GeoSpatialDistribution> geoSpatialDistributions = new ArrayList<GeoSpatialDistribution>();
        for (AbstractPolygon2D region : polygons) {
            geoSpatialDistributions.add(new GeoSpatialDistribution("test", region, false));
        }

        MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
                new MultiRegionGeoSpatialDistribution("multi", geoSpatialDistributions, false, overAllFill));

        Parameter.Default parameter = new Parameter.Default("coordinate", 2);
        likelihood.addData(parameter);
        double logL = likelihood.calculateLogLikelihood();

        System.out.println("\n(" + parameter.getParameterValue(0) + "," + parameter.getParameterValue(1) + ")");
        System.out.println("logL = " + logL);
        System.out.println("L = " + Math.exp(logL));

        parameter.setParameterValue(0, 11.40);
        parameter.setParameterValue(1, 106.90);

        likelihood.makeDirty();
        logL = likelihood.calculateLogLikelihood();

        System.out.println("\n(" + parameter.getParameterValue(0) + "," + parameter.getParameterValue(1) + ")");
        System.out.println("logL = " + logL);
        System.out.println("L = " + Math.exp(logL));


        //test reading in KML file with a single polygon without fillValue attribute
        polygons = Polygon2D.readKMLFile("Districts_polygons/21.035_106.063.kml");
        System.out.println("\n" + polygons.size() + " polygons found.");

        first = polygons.get(0);
        equalFills = first.hasFillValue();
        if (polygons.size() > 1) {
            for (int i = 1; i < polygons.size(); i++) {
                if (!(equalFills == polygons.get(i).hasFillValue())) {
                    System.out.println("Inconsistent fillValue attributes provided.");
                }
            }
        }

        if (equalFills) {
            double sum = 0.0;
            for (AbstractPolygon2D region : polygons) {
                sum += region.getFillValue();
            }
            if (Math.abs(sum - 1.0) > 1E-12) {
                System.out.println("\nFill values do not sum to 1 : " + sum);
            }
        }

        overAllFill = equalFills;
        geoSpatialDistributions = new ArrayList<GeoSpatialDistribution>();
        for (AbstractPolygon2D region : polygons) {
            geoSpatialDistributions.add(new GeoSpatialDistribution("test", region, false));
        }

        likelihood = new MultivariateDistributionLikelihood(
                new MultiRegionGeoSpatialDistribution("multi", geoSpatialDistributions, false, overAllFill));

        likelihood.addData(parameter);
        logL = likelihood.calculateLogLikelihood();

        System.out.println("\n(" + parameter.getParameterValue(0) + "," + parameter.getParameterValue(1) + ")");
        System.out.println("logL = " + logL);
        System.out.println("L = " + Math.exp(logL));

        parameter.setParameterValue(0, 21.0);
        parameter.setParameterValue(1, 106.03);

        likelihood.makeDirty();
        logL = likelihood.calculateLogLikelihood();

        System.out.println("\n(" + parameter.getParameterValue(0) + "," + parameter.getParameterValue(1) + ")");
        System.out.println("logL = " + logL);
        System.out.println("L = " + Math.exp(logL));

        parameter.setParameterValue(0, 21.0);
        parameter.setParameterValue(1, 106.10);

        likelihood.makeDirty();
        logL = likelihood.calculateLogLikelihood();

        System.out.println("\n(" + parameter.getParameterValue(0) + "," + parameter.getParameterValue(1) + ")");
        System.out.println("logL = " + logL);
        System.out.println("L = " + Math.exp(logL));

    }

}
