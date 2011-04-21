package dr.app.treespace;

import dr.app.treespace.CladeSystem.Clade;
import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreePlotter extends JComponent {

    public TreePlotter(final TreeLineages treeLineages) {
        setTrees(treeLineages);
    }

    public void setTrees(final TreeLineages treeLineages) {
        this.treeLineages = treeLineages;

        isCalibrated = false;
        repaint();
    }

    public void paint(Graphics g) {
//        if (trees == null || trees.size() == 0) {
//            return;
//        }
//
//        if (!isCalibrated) {
//            double[] yPosition = { 0.0 };
//
//            calibrateClade((Graphics2D)g, rootClade, yPosition);
//
//            scaleX = ((double)getWidth()) / rootClade.x;
//            scaleY = ((double)getHeight()) / yPosition[0];
//
//            isCalibrated = true;
//        }
//
//        paintClade((Graphics2D)g, rootClade);
    }

    private void paintClade(Graphics2D g2, Clade clade) {
        if (clade.children != null) {
            for (Clade child : clade.children.keySet()) {
                paintClade(g2, child);
                g2.setColor(Color.black);
                g2.setStroke(new BasicStroke(1.0F));
                g2.draw(new Line2D.Double(transformX(child.x), transformY(child.y), transformX(clade.x), transformY(clade.y)));
            }
        } else {
            paintTip(g2, clade);
        }

        float x = transformX(clade.x);
        float y = transformY(clade.y);

        Shape node = new Ellipse2D.Float(x, y, (float)(nodeWidth * scaleX), (float)(nodeHeight * scaleY));
        g2.setColor(Color.yellow);
        g2.fill(node);
        g2.setColor(Color.black);
        g2.draw(node);
    }

    private void paintTip(Graphics2D g2, Clade clade) {
        float x = transformX(clade.x);
        float y = transformY(clade.y);

        Shape node = new Ellipse2D.Float(x, y, (float)(nodeWidth * scaleX), (float)(nodeHeight * scaleY));
        g2.setColor(Color.yellow);
        g2.fill(node);
        g2.setColor(Color.black);
        g2.draw(node);
        if (clade.label != null) {
            g2.drawString(clade.label, x + xLabelOffset, yLabelOffset);
        }
    }

//    private double calibrate(Graphics2D g2, RootedTree tree, Node node, double[] yPosition) {
//        double y;
//        if (clade.children != null) {
//            y = 0.0;
//            for (Clade child : clade.children.keySet()) {
//                y += calibrateClade(g2, child, yPosition);
//            }
//
//            // todo: use an average weighted by support
//            y /= clade.children.keySet().size();
//        } else {
//            y = yPosition[0];
//            yPosition[0]++;
//        }
//        clade.x = (float)clade.bits.cardinality();
//        clade.y = (float)y;
//
//        return y;
//    }

    private float transformX(double x) {
        return (float)getWidth() - (float)(x * scaleX);
    }

    private float transformY(double y) {
        return (float)(y * scaleY);
    }

    private TreeLineages treeLineages;

    private double nodeWidth = 5.0;
    private double nodeHeight = 5.0;

    private float xLabelOffset = 5.0F;
    private float yLabelOffset = 5.0F;
    private boolean isCalibrated = false;

    private double scaleX = 1.0;
    private double scaleY = 1.0;

}
