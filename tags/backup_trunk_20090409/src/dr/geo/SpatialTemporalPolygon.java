package dr.geo;

import dr.xml.*;

/**
 * SpatialTemporalPolygon defines a spatial polygon in time.  This class could inherit from Polygon which uses ints
 * instead of reals to define space
 *
 * @author Marc A. Suchard
 * @author Alexei J. Drummond
 */


public class SpatialTemporalPolygon {

    public static final String POLYGON = "polygon";
    public static final String CLOSED = "closed";
    public static final String FILL_VALUE = "fillValue";

    public SpatialTemporalPolygon(double[] x, double[] y, double[] z, double[] t, boolean closed) {
        if (closed) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.t = t;
            length = x.length - 1;
        } else {
            this.x = wrapArray(x);
            this.y = wrapArray(y);
            this.z = wrapArray(z);
            this.t = wrapArray(t);
            length = x.length;
        }

        planarInZ = isPlanarInZ();

    }

    public void setFillValue(double value) {
        fillValue = value;
    }

    public double getFillValue() {
        return fillValue;
    }

    private boolean isPlanarInZ() {
        double z0 = z[0];
        boolean planar = true;
        for (int i = 0; i < length; i++) {
            if (z0 != z[i]) {
                planar = false;
                break;
            }
        }
        return planar;
    }

    private double[] wrapArray(double[] inArray) {
        if (inArray == null)
            return null;
        double[] outArray = new double[inArray.length + 1];
        System.arraycopy(inArray, 0, outArray, 0, inArray.length);
        outArray[inArray.length] = inArray[0];
        return outArray;
    }

    public boolean contains2DPoint(double inX, double inY) {

        if (!planarInZ)
            throw new RuntimeException("Only 2D polygons are currently implemented");
        boolean contains = false;

        // Take a horizontal ray from (inX,inY) to the right.
        // If ray across the polygon edges an odd # of times, the point is inside.
        for (int i = 0, j = length - 1; i < length; j = i++) {

             if ((((y[i] <= inY) && (inY < y[j])) ||
                    ((y[j] <= inY) && (inY < y[i]))) &&
                    (inX < (x[j] - x[i]) * (inY - y[i]) / (y[j] - y[i]) + x[i]))
                contains = !contains;
        }

        return contains;
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return POLYGON;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            KMLCoordinates coordinates = (KMLCoordinates) xo.getChild(KMLCoordinates.class);
            boolean closed = xo.getAttribute(CLOSED, false);

            if ((!closed && coordinates.length < 3) ||
                    (closed && coordinates.length < 4))
                throw new XMLParseException("Insufficient points in polygon '" + xo.getId() + "' to define a polygon in 2D");

            SpatialTemporalPolygon polygon = new SpatialTemporalPolygon(coordinates.x, coordinates.y, coordinates.z, null, closed);
            polygon.setFillValue(xo.getAttribute(FILL_VALUE,0.0));
            
            return polygon;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a polygon.";
        }

        public Class getReturnType() {
            return SpatialTemporalPolygon.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(KMLCoordinates.class),
                AttributeRule.newBooleanRule(CLOSED, true),
                AttributeRule.newDoubleRule(FILL_VALUE,true),
        };
    };

    protected double[] x;
    protected double[] y;
    protected double[] z;
    protected double[] t;

    protected int length;
    private boolean planarInZ;
    private double fillValue;

}
