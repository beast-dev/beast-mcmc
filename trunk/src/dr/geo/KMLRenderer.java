package dr.geo;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

/**
 * @author Alexei Drummond
 */
public class KMLRenderer {

    BufferedImage image;
    Rectangle2D bounds;

    java.util.List<Polygon2D> polygons;
    java.util.List<Shape> shapes;

    ViewTransform viewTransform;

    Color shapeColor;
    Color background;

    public KMLRenderer(String kmlFileName, Color shapeColor, Color background) {

        polygons = Polygon2D.readKMLFile(kmlFileName);

        this.shapeColor = shapeColor;
        this.background = background;

        System.out.println("Read " + polygons.size() + " polygons");

        System.out.println("Converting polygons to shapes");

        double minX = Double.MAX_VALUE, maxX = -Double.MAX_VALUE, minY = Double.MAX_VALUE, maxY = -Double.MAX_VALUE;

        shapes = new ArrayList<Shape>();

        for (Polygon2D p : polygons) {

            Shape shape = p.getShape();

            bounds = shape.getBounds();

            if (bounds.getMinX() < minX) minX = bounds.getMinX();
            if (bounds.getMaxX() > maxX) maxX = bounds.getMaxX();
            if (bounds.getMinY() < minY) minY = bounds.getMinY();
            if (bounds.getMaxY() > maxY) maxY = bounds.getMaxY();

            shapes.add(shape);
            System.out.print(".");
            System.out.flush();
        }

        bounds = new Rectangle2D.Double(minX, minY, maxX - minX, maxY - minY);

        System.out.println();

    }

    public Rectangle2D getBounds() {
        return bounds;
    }

    public void render(BufferedImage image) {
        Graphics2D g2d = image.createGraphics();

        g2d.setColor(background);
        g2d.fillRect(0, 0, image.getWidth(), image.getHeight());

        viewTransform = new ViewTransform(bounds, image.getWidth(), image.getHeight());

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2d.setColor(shapeColor);
        AffineTransform transform = viewTransform.getTransform();
        for (Shape s : shapes) {
            GeneralPath path = new GeneralPath(s);
            path.transform(transform);
            g2d.fill(path);
        }
    }
}
