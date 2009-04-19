package dr.geo;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateDistribution;
import dr.xml.*;

import java.util.List;
import java.util.ArrayList;

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

    public static final int dimPoint = 2; // Assumes 2D points only

    public GeoSpatialDistribution(SpatialTemporalPolygon region) {
        this.region = region;
    }

    public GeoSpatialDistribution(String label, SpatialTemporalPolygon region) {
        this.label = label;
        this.region = region;
    }

    public double logPdf(double[] x) {        
        if (region.contains2DPoint(x[0],x[1]))
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

    public String getLabel() { return label; }

    protected SpatialTemporalPolygon region;
    protected String label = null;


    public static XMLObjectParser FLAT_GEOSPATIAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FLAT_SPATIAL_DISTRIBUTION;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String label = xo.getAttribute(NODE_LABEL,"");

            List<GeoSpatialDistribution> geoSpatialDistributions = new ArrayList<GeoSpatialDistribution>();
            for(int i=0; i<xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof SpatialTemporalPolygon) {
                    SpatialTemporalPolygon region = (SpatialTemporalPolygon) xo.getChild(i);
                    geoSpatialDistributions.add(
                            new GeoSpatialDistribution(label,region)
                    );
                }
            }

            List<Parameter> parameters = new ArrayList<Parameter>();
            XMLObject cxo = (XMLObject) xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                Parameter spatialParameter = (Parameter) cxo.getChild(j);
                parameters.add(spatialParameter);
            }

            if (geoSpatialDistributions.size() == 1) {
                MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(geoSpatialDistributions.get(0));
                for(Parameter spatialParameter : parameters) {
                    if (spatialParameter.getDimension() != dimPoint)
                        throw new XMLParseException("Spatial priors currently only work in "+dimPoint+"D");
                    likelihood.addData(spatialParameter);
                }
                return likelihood;
            }

            if (parameters.size() == 1) {
                Parameter parameter = parameters.get(0);
                if (parameter.getDimension() % dimPoint != 0)
                    throw new XMLParseException("Spatial priors currently only work in "+dimPoint+"D");
                return new GeoSpatialCollectionModel(xo.getId(),parameter,geoSpatialDistributions);
            }

            throw new XMLParseException("Multiple separate parameters and multiple regions not yet implemented");

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(NODE_LABEL,true),
                new ElementRule(SpatialTemporalPolygon.class,1,Integer.MAX_VALUE),
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
