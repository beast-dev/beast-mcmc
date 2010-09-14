package dr.geo;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class MultiRegionGeoSpatialDistribution extends GeoSpatialDistribution {

    public MultiRegionGeoSpatialDistribution(String label) {
        super(label);
        regions = new ArrayList<GeoSpatialDistribution>();
        union = false;
    }

    public MultiRegionGeoSpatialDistribution(String label, List<GeoSpatialDistribution> regions, boolean union) {
        super(label);
        this.regions = regions;
        this.union = union;
    }

     public double logPdf(double[] x) {

        if (union) {
            for (GeoSpatialDistribution region : regions) {
                if (region.logPdf(x) == 0.0) { // matches
                    return 0.0;
                }
            }
            return Double.NEGATIVE_INFINITY;
        } // else is intersection

        for (GeoSpatialDistribution region : regions) {
            if (region.logPdf(x) == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0.0;
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

    public boolean getUnion() {
        return union;
    }

    private final List<GeoSpatialDistribution> regions;
    private final boolean union;

}
