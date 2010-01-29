package dr.geo;

import dr.math.distributions.MultivariateNormalDistribution;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class VisualizeKMLBrownianBridge extends VisualizeBrownianBridge2D {

    List<Polygon2D> polygons;

    public VisualizeKMLBrownianBridge(String kmlFileName) {

        polygons = Polygon2D.readKMLFile(kmlFileName);

        System.out.println("Read " + polygons.size() + " polygons");

        start = new SpaceTime(0, new double[]{0, 45});
        end = new SpaceTime(1, new double[]{21, 40});

        topLeft = new Point2D.Double(-5, 8);
        bottomRight = new Point2D.Double(50, 70);


        System.out.println("Converting polygons to shapes");
        shapes = new ArrayList<Shape>();
        for (Polygon2D p : polygons) {
            shapes.add(getShape(p));
            System.out.print(".");
            System.out.flush();
        }
        System.out.println();

        rejector = new SpaceTimeRejector() {

            public boolean reject(SpaceTime point) {
                Point2D p = new Point2D.Double(point.getX(0), point.getX(1));
                for (Shape s : shapes) {
                    if (s.contains(p)) return true;
                }
                return false;
            }
        };

        mnd = new MultivariateNormalDistribution(new double[]{0.0}, new double[][]{{0.1, 0}, {0, 0.1}});

        shapeColor = Color.BLUE;

        if (rejector.reject(start) || rejector.reject(end)) {
            throw new RuntimeException("Start or end in water");
        }
    }

    Shape getShape(Polygon2D poly) {
        GeneralPath path = new GeneralPath();

        LinkedList<Point2D> points = poly.point2Ds;
        path.moveTo((float) points.get(0).getX(), (float) points.get(0).getY());

        System.out.println("x=" + points.get(0).getX() + ", y=" + points.get(0).getY());

        for (int i = 1; i < points.size(); i++) {
            path.lineTo((float) points.get(i).getX(), (float) points.get(i).getY());
        }
        path.closePath();
        return path;
    }

    AffineTransform getTranslate() {
        return AffineTransform.getTranslateInstance(-topLeft.getX(), -bottomRight.getY());
    }

    AffineTransform getScale() {
        return AffineTransform.getScaleInstance(scaleX, -scaleY);
    }


    public static void main(String[] args) {

        JFrame frame = new JFrame("Europe");
        frame.getContentPane().add(BorderLayout.CENTER, new VisualizeKMLBrownianBridge(args[0]));
        frame.setSize(600, 600);
        frame.setVisible(true);
    }

}
