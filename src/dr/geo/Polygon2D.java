package dr.geo;

import dr.xml.*;

import java.awt.geom.Rectangle2D;
import java.awt.geom.Point2D;
import java.util.LinkedList;


/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class Polygon2D {
    public static final String POLYGON = "polygon";
    public static final String CLOSED = "closed";
    public static final String FILL_VALUE = "fillValue";

    public Polygon2D(LinkedList<Point2D> points, boolean closed) {
        this.Point2Ds = points;
        if (!closed) {
            Point2D start = points.get(0);
            points.add(start);
        }
        length = points.size() - 1;
    }

    public Polygon2D() {
        length = 0;
        Point2Ds = new LinkedList<Point2D>();
    }

    public void addPoint2D(Point2D Point2D) {
        if (Point2Ds.size() == 0)
            Point2Ds.add(Point2D);
        else if (Point2Ds.size() == 1) {
            Point2Ds.add(Point2D);
            Point2Ds.add(Point2Ds.get(0));
        } else {
            Point2D last = Point2Ds.removeLast();
            Point2Ds.add(Point2D);
            if(!last.equals(Point2D))
                Point2Ds.add(last);
        }

        length = Point2Ds.size() - 1;
    }

    public void setFillValue(double value) {
        fillValue = value;
    }

    public double getFillValue() {
        return fillValue;
    }

    public boolean containsPoint2D(Point2D Point2D) {

        final double inX = Point2D.getX();
        final double inY = Point2D.getY();
        boolean contains = false;

        // Take a horizontal ray from (inX,inY) to the right.
        // If ray across the polygon edges an odd # of times, the Point2D is inside.
        for (int i = 0, j = length - 1; i < length; j = i++) {

            final double xi = Point2Ds.get(i).getX();
            final double xj = Point2Ds.get(j).getX();

            final double yi = Point2Ds.get(i).getY();
            final double yj = Point2Ds.get(j).getY();

            if ((((yi <= inY) && (inY < yj)) ||
                    ((yj <= inY) && (inY < yi))) &&
                    (inX < (xj - xi) * (inY - yi) / (yj - yi) + xi))
                contains = !contains;
        }

        return contains;
    }

    private enum Side {
        left, right, top, bottom
    }

    public Polygon2D clip(Rectangle2D boundingBox) {

        LinkedList<Point2D> clippedPolygon = new LinkedList<Point2D>();

        Point2D p;        // current Point2D
        Point2D p2;       // next Point2D

        // make copy of original polygon to work with
        LinkedList<Point2D> workPoly = new LinkedList<Point2D>(Point2Ds);

        // loop through all for clipping edges
        for (Side side : Side.values()) {
            clippedPolygon.clear();
            for (int i = 0; i < workPoly.size() - 1; i++) {
                p = workPoly.get(i);
                p2 = workPoly.get(i + 1);
                if (isInsideClip(p, side, boundingBox)) {
                    if (isInsideClip(p2, side, boundingBox))
                        // here both Point2Ds are inside the clipping window so add the second one
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
        return new Polygon2D(clippedPolygon,true);
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
        sb.append("polygon[\n");
        for(Point2D pt : Point2Ds) {
            sb.append("\t");
            sb.append(pt);
            sb.append("\n");
        }
        sb.append("]");
        return sb.toString();
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
                throw new XMLParseException("Insufficient Point2Ds in polygon '" + xo.getId() + "' to define a polygon in 2D");

            LinkedList<Point2D> Point2Ds = new LinkedList<Point2D>();
            for (int i = 0; i < coordinates.length; i++)
                Point2Ds.add(new Point2D.Double(coordinates.x[i], coordinates.y[i]));

            Polygon2D polygon = new Polygon2D(Point2Ds, closed);
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
            return SpatialTemporalPolygon.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(KMLCoordinates.class),
                AttributeRule.newBooleanRule(CLOSED, true),
                AttributeRule.newDoubleRule(FILL_VALUE, true),
        };
    };

    public static void main(String[] args) {
        Polygon2D polygon = new Polygon2D();
        polygon.addPoint2D(new Point2D.Double(-10,-10));
        polygon.addPoint2D(new Point2D.Double(-10,50));
        polygon.addPoint2D(new Point2D.Double(10,50));
        polygon.addPoint2D(new Point2D.Double(10,-10));
        System.out.println(polygon);
        System.out.println("");

        Point2D pt = new Point2D.Double(0,0);
        System.out.println("polygon contains "+pt+": "+polygon.containsPoint2D(pt));
        pt = new Point2D.Double(100,100);
        System.out.println("polygon contains "+pt+": "+polygon.containsPoint2D(pt));
        System.out.println("");

        Rectangle2D boundingBox = new Rectangle2D.Double(0,0,100,100);  // defines lower-left corner and width/height
        System.out.println(boundingBox);
        Polygon2D myClip = polygon.clip(boundingBox);
        System.out.println(myClip);

    }

    protected LinkedList<Point2D> Point2Ds;
    protected int length;
    private double fillValue;

}