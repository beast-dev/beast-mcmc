package dr.geo;

import dr.math.distributions.MultivariateNormalDistribution;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class VisualizeBrownianBridge2D extends JComponent {

    MultivariateNormalDistribution mnd;
    SpaceTime start, end;
    SpaceTimeRejector rejector;

    List<Shape> shapes;

    Paint shapeColor = Color.GRAY;

    Point2D topLeft;
    Point2D bottomRight;

    double scaleX;
    double scaleY;


    public VisualizeBrownianBridge2D() {

        start = new SpaceTime(0, new double[]{0, 0});
        end = new SpaceTime(1, new double[]{1, 1});

        topLeft = new Point2D.Double(-0.2, -0.2);
        bottomRight = new Point2D.Double(1.2, 1.2);

        shapes = new ArrayList<Shape>();
        shapes.add(new Ellipse2D.Double(0.25, 0.25, 0.4, 0.4));
        shapes.add(new Ellipse2D.Double(0.5, 0.7, 0.2, 0.2));
        shapes.add(new Ellipse2D.Double(0.8, 0.2, 0.15, 0.15));

        rejector = new SpaceTimeRejector() {

            public boolean reject(SpaceTime point) {
                Point2D p = new Point2D.Double(point.getX(0), point.getX(1));
                for (Shape s : shapes) {
                    if (s.contains(p)) return true;
                }
                return false;
            }
        };

        mnd = new MultivariateNormalDistribution(new double[]{0.0}, new double[][]{{6, 0}, {0, 6}});
    }

    public void setShapeColor(Color c) {
        shapeColor = c;
    }

    void computeScales() {
        scaleX = getWidth() / (bottomRight.getX() - topLeft.getX());
        scaleY = getHeight() / (bottomRight.getY() - topLeft.getY());
    }

    public void paintComponent(Graphics g) {

        System.out.println("entering paintComponent()");

        computeScales();

        Graphics2D g2d = (Graphics2D) g;

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2d.setStroke(new BasicStroke(1.5f));

        System.out.println("Painting shapes");
        for (Shape s : shapes) {
            System.out.print(".");
            System.out.flush();
            GeneralPath path = new GeneralPath(s);
            path.transform(getFullTransform());
            g2d.setPaint(shapeColor);
            g2d.fill(path);
        }


        AffineTransform transform = getFullTransform();

        for (int r = 0; r < 1; r++) {
            g2d.setPaint(new Color((float) Math.random(), (float) Math.random(), (float) Math.random()));

            List<SpaceTime> points = MultivariateBrownianBridge.divideConquerBrownianBridge(mnd, start, end, 6, rejector);

            GeneralPath path = new GeneralPath();
            path.moveTo(points.get(0).getX(0), points.get(0).getX(1));
            //System.out.println(points.get(0));
            for (int i = 1; i < points.size(); i++) {
                path.lineTo(points.get(i).getX(0), points.get(i).getX(1));
                //System.out.println(points.get(i));
            }

            path.transform(getFullTransform());

            g2d.draw(path);
        }

        g2d.setPaint(Color.black);

        Point2D startPointRaw = new Point2D.Double(start.getX(0), start.getX(1));
        Point2D endPointRaw = new Point2D.Double(end.getX(0), end.getX(1));
        Point2D startPointT = new Point2D.Double();
        Point2D endPointT = new Point2D.Double();

        transform.transform(startPointRaw, startPointT);
        transform.transform(endPointRaw, endPointT);

        Shape startShape = new Ellipse2D.Double(startPointT.getX() - 3, startPointT.getY() - 3, 6, 6);
        Shape endShape = new Ellipse2D.Double(endPointT.getX() - 3, endPointT.getY() - 3, 6, 6);

        g2d.fill(startShape);
        g2d.fill(endShape);


        System.out.println("leaving paintComponent()");

    }

    AffineTransform getScale() {
        return AffineTransform.getScaleInstance(scaleX, scaleY);
    }

    AffineTransform getTranslate() {
        return AffineTransform.getTranslateInstance(-topLeft.getX(), -topLeft.getY());
    }

    AffineTransform getFullTransform() {
        AffineTransform transform = getScale();
        transform.concatenate(getTranslate());
        return transform;
    }

    public static void main(String[] args) {

        JFrame frame = new JFrame("Boulders");
        frame.getContentPane().add(BorderLayout.CENTER, new VisualizeBrownianBridge2D());
        frame.setSize(600, 600);
        frame.setVisible(true);
    }
}
