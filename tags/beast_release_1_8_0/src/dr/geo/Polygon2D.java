package dr.geo;

import dr.geo.cartogram.CartogramMapping;
import dr.xml.*;
import dr.util.HeapSort;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;


/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class Polygon2D {
    public static final String POLYGON = "polygon";
    public static final String CLOSED = "closed";
    public static final String FILL_VALUE = "fillValue";
    public static final String CIRCLE = "circle";
    public static final String NUMBER_OF_POINTS = "numberOfPoints";
    public static final String RADIUS = "radius";
    public static final String CENTER = "center";
    public static final String LATITUDE = "latitude";
    public static final String LONGITUDE = "longitude";

    public Polygon2D(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new RuntimeException("Unbalanced arrays");
        }

        if (x[0] != x[x.length - 1] && y[0] != y[y.length - 1]) {
            double[] newX = new double[x.length + 1];
            double[] newY = new double[y.length + 1];
            System.arraycopy(x, 0, newX, 0, x.length);
            System.arraycopy(y, 0, newY, 0, y.length);
            newX[x.length] = x[0];
            newY[y.length] = y[0];
            this.x = newX;
            this.y = newY;
        } else {
            this.x = x;
            this.y = y;
        }
        length = this.x.length - 1;

    }

    public Polygon2D(List<Point2D> points, boolean closed) {
        this.point2Ds = points;
        if (!closed) {
            Point2D start = points.get(0);
            points.add(start);
        }
        convertPointsToArrays();
        length = points.size() - 1;
    }

    public Polygon2D() {
        length = 0;
        point2Ds = new ArrayList<Point2D>();
    }

    public String getID() {
        return id;
    }

    public Polygon2D(Element e) {

//        System.err.println("parsing polygon");

        List<Element> children = e.getChildren();
        id = e.getAttributeValue(XMLParser.ID);

        parseCoordinates(e);
    }

    private void parseCoordinates(Element element) {

        if (element.getName().equalsIgnoreCase(KMLCoordinates.COORDINATES)) {
            String value = element.getTextTrim();
            StringTokenizer st1 = new StringTokenizer(value, KMLCoordinates.POINT_SEPARATORS);
            int count = st1.countTokens();
//            System.out.println(count + " tokens");

            point2Ds = new ArrayList<Point2D>(count);
            for (int i = 0; i < count; i++) {
                String line = st1.nextToken();
                StringTokenizer st2 = new StringTokenizer(line, KMLCoordinates.SEPARATOR);
                if (st2.countTokens() < 2 || st2.countTokens() > 3)
                    throw new IllegalArgumentException("All KML coordinates must contain (X,Y) or (X,Y,Z) values.  Error in element '" + line + "'");
                final double x = Double.valueOf(st2.nextToken());
                final double y = Double.valueOf(st2.nextToken());
                point2Ds.add(new Point2D.Double(x, y));
            }
            convertPointsToArrays();
            length = point2Ds.size() - 1;
        } else {
            for (Object child : element.getChildren()) {

                if (child instanceof Element) {
                    parseCoordinates((Element) child);
                }
            }
        }
    }

    Shape getShape() {
        GeneralPath path = new GeneralPath();

        List<Point2D> points = point2Ds;
        path.moveTo((float) points.get(0).getX(), (float) points.get(0).getY());

        for (int i = 1; i < points.size(); i++) {
            path.lineTo((float) points.get(i).getX(), (float) points.get(i).getY());
        }
        path.closePath();
        return path;
    }

    private void convertPointsToArrays() {
        final int length = point2Ds.size();
        if (x == null || x.length != length) {
            x = new double[length];
            y = new double[length];
        }
        Iterator<Point2D> it = point2Ds.iterator();
        for (int i = 0; i < length; i++) {
            final Point2D point = it.next();
            x[i] = point.getX();
            y[i] = point.getY();
        }
    }

    public void addPoint2D(Point2D point2D) {
        if (point2Ds.size() == 0)
            point2Ds.add(point2D);
        else if (point2Ds.size() == 1) {
            point2Ds.add(point2D);
            point2Ds.add(point2Ds.get(0));
        } else {
            Point2D last = point2Ds.remove(point2Ds.size() - 1);
            point2Ds.add(point2D);
            if (!last.equals(point2D))
                point2Ds.add(last);
        }
        convertPointsToArrays();
        length = point2Ds.size() - 1;
    }

    public Point2D getPoint2D(int x) {
        if (x > length +1) {
            throw new RuntimeException("Polygon only has length"+length);
        } else {
            return point2Ds.get(x);
        }
    }

    public boolean containsPoint2D(Point2D Point2D) {

        final double inX = Point2D.getX();
        final double inY = Point2D.getY();
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

    public boolean bordersPoint2D(Point2D Point2D) {
        boolean borders = false;

        Iterator<Point2D> it = point2Ds.iterator();
        for (int i = 0; i < length; i++) {
            Point2D point = it.next();
            if (point.equals(Point2D)) {
                borders = true;
            }
        }
        return borders;
    }

    public void setFillValue(double value) {
        fillValue = value;
    }

    public double getFillValue() {
        return fillValue;
    }

    public double getLength() {
        return length;
    }

//    public boolean containsPoint2D(Point2D Point2D) { // this takes 3 times as long as the above code, why???
//
//        final double inX = Point2D.getX();
//        final double inY = Point2D.getY();
//        boolean contains = false;
//
//        // Take a horizontal ray from (inX,inY) to the right.
//        // If ray across the polygon edges an odd # of times, the Point2D is inside.
//
//        final Point2D end   = point2Ds.get(length-1); // assumes closed
//        double xi = end.getX();
//        double yi = end.getY();
//
//        Iterator<Point2D> listIterator = point2Ds.iterator();
//
//        for(int i=0; i<length; i++) {
//
//            final double xj = xi;
//            final double yj = yi;
//
//            final Point2D next = listIterator.next();
//            xi = next.getX();
//            yi = next.getY();
//
//            if ((((yi <= inY) && (inY < yj)) ||
//                    ((yj <= inY) && (inY < yi))) &&
//                    (inX < (xj - xi) * (inY - yi) / (yj - yi) + xi))
//                contains = !contains;
//        }
//        return contains;
//    }

    private enum Side {
        left, right, top, bottom
    }

    public Polygon2D clip(Rectangle2D boundingBox) {

        LinkedList<Point2D> clippedPolygon = new LinkedList<Point2D>();

        Point2D p;        // current Point2D
        Point2D p2;       // next Point2D

        // make copy of original polygon to work with
        LinkedList<Point2D> workPoly = new LinkedList<Point2D>(point2Ds);

        // loop through all for clipping edges
        for (Side side : Side.values()) {
            clippedPolygon.clear();
            for (int i = 0; i < workPoly.size() - 1; i++) {
                p = workPoly.get(i);
                p2 = workPoly.get(i + 1);
                if (isInsideClip(p, side, boundingBox)) {
                    if (isInsideClip(p2, side, boundingBox))
                        // here both point2Ds are inside the clipping window so add the second one
                        clippedPolygon.add(p2);
                    else
                        // the seond Point2D is outside so add the intersection Point2D
                        clippedPolygon.add(intersectionPoint2D(side, p, p2, boundingBox));
                } else {
                    // so first Point2D is outside the window here
                    if (isInsideClip(p2, side, boundingBox)) {
                        // the following Point2D is inside so add the insection Point2D and also p2
                        clippedPolygon.add(intersectionPoint2D(side, p, p2, boundingBox));
                        clippedPolygon.add(p2);
                    }
                }
            }
            // make sure that first and last element are the same, we want a closed polygon
            if (!clippedPolygon.getFirst().equals(clippedPolygon.getLast()))
                clippedPolygon.add(clippedPolygon.getFirst());
            // we have to keep on working with our new clipped polygon
            workPoly = new LinkedList<Point2D>(clippedPolygon);
        }
        return new Polygon2D(clippedPolygon, true);
    }

    public void transformByMapping(CartogramMapping mapping) {
        for (int i = 0; i < length + 1; i++) {
            point2Ds.set(i, mapping.map(point2Ds.get(i)));
        }
        convertPointsToArrays();
    }

    public void swapXYs() {
        for (int i = 0; i < length + 1; i++) {
            point2Ds.set(i, new Point2D.Double(point2Ds.get(i).getY(), point2Ds.get(i).getX()));
        }
        convertPointsToArrays();
    }

    public void rescale(double longMin, double longwidth, double gridXSize, double latMax, double latwidth, double gridYSize) {
        for (int i = 0; i < length + 1; i++) {
            point2Ds.set(i, new Point2D.Double(((point2Ds.get(i).getX() - longMin) * (gridXSize / longwidth)), ((latMax - point2Ds.get(i).getY()) * (gridYSize / latwidth))));
        }
        convertPointsToArrays();
    }

    public void rescaleToPositiveCoordinates() {

        double[][] xyMinMax = getXYMinMax();
        double shiftX = 0;
        double shiftY = 0;

        if (xyMinMax[0][0] < 0){
            shiftX = -xyMinMax[0][0];
        }
        if (xyMinMax[1][0] < 0){
            shiftY = -xyMinMax[1][0];
        }

        if ((shiftX < 0) || (shiftY < 0)) {
            for (int i = 0; i < length + 1; i++) {
                point2Ds.set(i, new Point2D.Double(point2Ds.get(i).getX()+shiftX, point2Ds.get(i).getY()+shiftY));
            }
            convertPointsToArrays();
        }
    }

    public double[][] getXYMinMax(){

        int[] indicesX = new int[x.length];
        int[] indicesY = new int[y.length];
        HeapSort.sort(x, indicesX);
        HeapSort.sort(y, indicesY);

        double[][] returnArray = new double[2][2];
        returnArray[0][0] = x[indicesX[0]];
        returnArray[0][1] = x[indicesX[indicesX.length - 1]];
        returnArray[1][0] = y[indicesY[0]];
        returnArray[1][1] = y[indicesY[indicesY.length - 1]];

        return returnArray;
    }


    // Here is a formula for the area of a polygon with vertices {(xk,yk): k = 1,...,n}:
    //   Area = 1/2 [(x1*y2 - x2*y1) + (x2*y3 - x3*y2) + ... + (xn*y1 - x1*yn)].
    //   This formula appears in an Article by Gil Strang of MIT
    //   on p. 253 of the March 1993 issue of The American Mathematical Monthly, with the note that it is
    //   "known, but not well known". There is also a very brief discussion of proofs and other references,
    //   including an article by Bart Braden of Northern Kentucky U., a known Mathematica enthusiast.
    public double calculateArea() {

//        rescaleToPositiveCoordinates();

        double area = 0;
        //we can implement it like this because the polygon is closed (point2D.get(0) = point2D.get(length + 1)
        for (int i = 0; i < length; i++) {
            area += (x[i] * y[i + 1] - x[i + 1] * y[i]);
        }

        return (Math.abs(area / 2));
    }

    public Point2D getCentroid() {

//        rescaleToPositiveCoordinates();

        Point2D centroid = new Point2D.Double();
        double area = calculateArea();
        double cx=0,cy=0;
        
        double factor;

        //we can implement it like this because the polygon is closed (point2D.get(0) = point2D.get(length + 1)
        for (int i = 0; i < length; i++) {
            factor = (x[i] * y[i + 1] - x[i + 1] * y[i]);
            cx += (x[i] * x[i + 1])*factor;
            cy += (y[i] * y[i + 1])*factor;
        }
        double constant = 1/(area*6);
        cx*=constant;
        cy*=constant;
        centroid.setLocation(cx,cy);
        System.out.println("centroid = "+cx+","+cy);
        return centroid;
    }

    private static LinkedList<Point2D> getCirclePoints(double centerLat, double centerLong, int numberOfPoints, double radius) {

        LinkedList<Point2D> Point2Ds = new LinkedList<Point2D>();

        double lat1, long1;
        double d_rad;
        double delta_pts;
        double radial, lat_rad, dlon_rad, lon_rad;

        // convert coordinates to radians
        lat1 = Math.toRadians(centerLat);
        long1 = Math.toRadians(centerLong);

        //radius is in meters
        d_rad = radius / 6378137;

        // loop through the array and write points
        for (int i = 0; i <= numberOfPoints; i++) {
            delta_pts = 360 / (double) numberOfPoints;
            radial = Math.toRadians((double) i * delta_pts);

            //This algorithm is limited to distances such that dlon < pi/2
            lat_rad = Math.asin(Math.sin(lat1) * Math.cos(d_rad) + Math.cos(lat1) * Math.sin(d_rad) * Math.cos(radial));
            dlon_rad = Math.atan2(Math.sin(radial) * Math.sin(d_rad) * Math.cos(lat1), Math.cos(d_rad) - Math.sin(lat1) * Math.sin(lat_rad));
            lon_rad = ((long1 + dlon_rad + Math.PI) % (2 * Math.PI)) - Math.PI;

            Point2Ds.add(new Point2D.Double(Math.toDegrees(lat_rad), Math.toDegrees(lon_rad)));

        }
        return Point2Ds;
    }

    private static boolean isInsideClip(Point2D p, Side side, Rectangle2D boundingBox) {
        if (side == Side.top)
            return (p.getY() <= boundingBox.getMaxY());
        else if (side == Side.bottom)
            return (p.getY() >= boundingBox.getMinY());
        else if (side == Side.left)
            return (p.getX() >= boundingBox.getMinX());
        else if (side == Side.right)
            return (p.getX() <= boundingBox.getMaxX());
        else
            throw new RuntimeException("Error in Polygon");
    }


    private static Point2D intersectionPoint2D(Side side, Point2D p1, Point2D p2, Rectangle2D boundingBox) {

        if (side == Side.top) {
            final double topEdge = boundingBox.getMaxY();
            return new Point2D.Double(p1.getX() + (topEdge - p1.getY()) * (p2.getX() - p1.getX()) / (p2.getY() - p1.getY()), topEdge);
        } else if (side == Side.bottom) {
            final double bottomEdge = boundingBox.getMinY();
            return new Point2D.Double(p1.getX() + (bottomEdge - p1.getY()) * (p2.getX() - p1.getX()) / (p2.getY() - p1.getY()), bottomEdge);
        } else if (side == Side.right) {
            final double rightEdge = boundingBox.getMaxX();
            return new Point2D.Double(rightEdge, p1.getY() + (rightEdge - p1.getX()) * (p2.getY() - p1.getY()) / (p2.getX() - p1.getX()));
        } else if (side == Side.left) {
            final double leftEdge = boundingBox.getMinX();
            return new Point2D.Double(leftEdge, p1.getY() + (leftEdge - p1.getX()) * (p2.getY() - p1.getY()) / (p2.getX() - p1.getX()));
        }
        return null;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(POLYGON).append("[\n");
        for (Point2D pt : point2Ds) {
            sb.append("\t");
            sb.append(pt);
            sb.append("\n");
        }
        sb.append("]");
        return sb.toString();
    }

    public static void readKMLElement(Element element, List<Polygon2D> polygons) {

        if (element.getName().equalsIgnoreCase(POLYGON)) {
            Polygon2D polygon = new Polygon2D(element);
            polygons.add(polygon);
        } else {
            for (Object child : element.getChildren()) {

                if (child instanceof Element) {
                    readKMLElement((Element) child, polygons);
                }
            }

        }
    }

    public static List<Polygon2D> readKMLFile(String fileName) {

        List<Polygon2D> polygons = new ArrayList<Polygon2D>();
        try {

            SAXBuilder builder = new SAXBuilder();
            builder.setValidation(false);
            builder.setIgnoringElementContentWhitespace(true);
            Document doc = builder.build(new File(fileName));
            Element root = doc.getRootElement();
            if (!root.getName().equalsIgnoreCase("KML"))
                throw new RuntimeException("Not a KML file");

            readKMLElement(root, polygons);

        } catch (IOException e) {
            e.printStackTrace();
        } catch (JDOMException e) {
            e.printStackTrace();
        }
        return polygons;
    }
//
//    public Element toXML() {
//        return new KMLCoordinates(x,y).toXML();
//    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return POLYGON;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            LinkedList<Point2D> Point2Ds = new LinkedList<Point2D>();
            boolean closed;
            Polygon2D polygon;

            if (xo.getChild(Polygon2D.class) != null) { // This is a regular polygon

                polygon = (Polygon2D) xo.getChild(Polygon2D.class);

            } else { // This is an arbitrary polygon

                KMLCoordinates coordinates = (KMLCoordinates) xo.getChild(KMLCoordinates.class);
                closed = xo.getAttribute(CLOSED, false);

                if ((!closed && coordinates.length < 3) ||
                        (closed && coordinates.length < 4))
                    throw new XMLParseException("Insufficient point2Ds in polygon '" + xo.getId() + "' to define a polygon in 2D");

                for (int i = 0; i < coordinates.length; i++)
                    Point2Ds.add(new Point2D.Double(coordinates.x[i], coordinates.y[i]));

                polygon = new Polygon2D(Point2Ds, closed);

            }

            polygon.setFillValue(xo.getAttribute(FILL_VALUE, 0.0));

            return polygon;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a polygon.";
        }

        public Class getReturnType() {
            return Polygon2D.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new XORRule(
                        new ElementRule(KMLCoordinates.class),
                        new ElementRule(Polygon2D.class)
                ),
                AttributeRule.newBooleanRule(CLOSED, true),
                AttributeRule.newDoubleRule(FILL_VALUE, true),
        };
    };

    public static XMLObjectParser CIRCLE_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CIRCLE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double latitude = xo.getDoubleAttribute(LATITUDE);
            double longitude = xo.getDoubleAttribute(LONGITUDE);
            double radius = xo.getDoubleAttribute(RADIUS);
            int num = xo.getAttribute(NUMBER_OF_POINTS, 50); // default = 50

            LinkedList<Point2D> Point2Ds = getCirclePoints(latitude,
                    longitude,
                    num,
                    radius);

            return new Polygon2D(Point2Ds, true);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a regular circle polygon.";
        }

        public Class getReturnType() {
            return Polygon2D.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(LATITUDE),
                AttributeRule.newDoubleRule(LONGITUDE),
                AttributeRule.newDoubleRule(RADIUS),
                AttributeRule.newIntegerRule(NUMBER_OF_POINTS, true),
        };
    };

    public static void main(String[] args) {
        Polygon2D polygon = new Polygon2D();
        polygon.addPoint2D(new Point2D.Double(-10, -10));
        polygon.addPoint2D(new Point2D.Double(-10, 50));
        polygon.addPoint2D(new Point2D.Double(10, 50));
        polygon.addPoint2D(new Point2D.Double(10, -10));
        System.out.println(polygon);
        System.out.println("");

        Point2D pt = new Point2D.Double(0, 0);
        System.out.println("polygon contains " + pt + ": " + polygon.containsPoint2D(pt));
        pt = new Point2D.Double(100, 100);
        System.out.println("polygon contains " + pt + ": " + polygon.containsPoint2D(pt));
        System.out.println("");

        Rectangle2D boundingBox = new Rectangle2D.Double(0, 0, 100, 100);  // defines lower-left corner and width/height
        System.out.println(boundingBox);
        Polygon2D myClip = polygon.clip(boundingBox);
        System.out.println(myClip);

    }

    protected List<Point2D> point2Ds;

    protected int length;
    private double fillValue;
    private String id;

    protected double[] x;
    protected double[] y;


}