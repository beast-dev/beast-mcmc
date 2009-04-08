package dr.geo;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateDistribution;
import dr.xml.*;

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

            SpatialTemporalPolygon region = (SpatialTemporalPolygon) xo.getChild(SpatialTemporalPolygon.class);

            MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(new GeoSpatialDistribution(label,region));

            XMLObject cxo = (XMLObject) xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                Parameter spatialParameter = (Parameter) cxo.getChild(j);
                if (spatialParameter.getDimension() != 2)
                    throw new RuntimeException("Spatial priors currently only work in 2D");
                likelihood.addData(spatialParameter);
            }
            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(NODE_LABEL,true),
                new ElementRule(SpatialTemporalPolygon.class),
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
