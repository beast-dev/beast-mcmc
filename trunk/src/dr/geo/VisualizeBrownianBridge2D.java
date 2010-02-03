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

            ArrayList<Reject> rejects = new ArrayList<Reject>();

            public boolean reject(SpaceTime point, int attribute) {
                Point2D p = new Point2D.Double(point.getX(0), point.getX(1));
                for (Shape s : shapes) {
                    if (s.contains(p)) {
                        rejects.add(new Reject(attribute, point));
                        return true;
                    }
                }
                return false;
            }

            // removes all rejects
            public void reset() {
                rejects.clear();
            }

            public List<Reject> getRejects() {
                return rejects;
            }

//            private boolean stop = false;
//
//            public boolean getStop() {
//                return stop;
//            }
//
//            public void setStop(boolean stop) {
//                this.stop = stop;
//            }


        };

        mnd = new MultivariateNormalDistribution(new double[]{0.0}, new double[][]{{10, 0}, {0, 10}});
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

        for (int r = 0; r < getTrials(); r++) {
            g2d.setPaint(new Color((float) Math.random(), (float) Math.random(), (float) Math.random()));

            List<SpaceTime> points = null;

            int topLevelRejects = -1;
            while (points == null) {
                topLevelRejects += 1;
                points = MultivariateBrownianBridge.divideConquerBrownianBridge(mnd, start, end, getMaxDepth(), getMaxTries(), rejector);
            }

            Paint old = g2d.getPaint();
            g2d.setPaint(Color.yellow);

            String rejectString = computeRejectString(rejector, topLevelRejects);
            g2d.drawString(rejectString, 10, getHeight() - 20);
            //rejector.reset();
            g2d.setPaint(old);

            GeneralPath path = new GeneralPath();
            path.moveTo((float) points.get(0).getX(0), (float) points.get(0).getX(1));
            //System.out.println(points.get(0));
            for (int i = 1; i < points.size(); i++) {
                path.lineTo((float) points.get(i).getX(0), (float) points.get(i).getX(1));
                //System.out.println(points.get(i));
            }

            path.transform(getFullTransform());

            g2d.draw(path);
        }

        g2d.setPaint(Color.black);

        paintDot(start, 3, transform, g2d);
        paintDot(end, 3, transform, g2d);

        System.out.println("leaving paintComponent()");

    }

    private String computeRejectString(SpaceTimeRejector rejector, int topLevelRejects) {

        int[] rejectCounts = new int[9];
        for (Reject r : rejector.getRejects()) {
            rejectCounts[Math.min(r.getDepth() - 1, rejectCounts.length - 1)] += 1;
        }
        StringBuilder builder = new StringBuilder();
        builder.append("Rejects top-level=" + topLevelRejects);
        for (int i = 0; i < rejectCounts.length; i++) {
            builder.append("  " + (i + 1) + ":" + rejectCounts[i]);
        }
        return builder.toString();
    }

    void paintDot(SpaceTime s, double radius, AffineTransform transform, Graphics2D g2d) {

        Point2D pointRaw = new Point2D.Double(s.getX(0), s.getX(1));
        Point2D pointT = new Point2D.Double();

        transform.transform(pointRaw, pointT);

        Shape pointShape = new Ellipse2D.Double(pointT.getX() - radius, pointT.getY() - radius, 2.0 * radius, 2.0 * radius);

        g2d.fill(pointShape);

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

    public int getMaxDepth() {
        return 10;
    }

    public int getTrials() {

        return 10;
    }

    public int getMaxTries() {
        return 10;
    }
}
