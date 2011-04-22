package dr.app.treespace;

import dr.app.treespace.CladeSystem.Clade;
import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;
import jebl.math.Random;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.image.ColorModel;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreePlotter extends JComponent {

    public TreePlotter() {
        setTrees(treeLineages);
    }

    public TreePlotter(final TreeLineages treeLineages) {
        setTrees(treeLineages);
    }

    public void setTrees(final TreeLineages treeLineages) {
        this.treeLineages = treeLineages;

        isCalibrated = false;
        repaint();
    }

    public void paint(Graphics g) {
        if (treeLineages == null || treeLineages.getRootLineages().size() == 0) {
            return;
        }

        if (!isCalibrated) {
            scaleX = getBounds().getWidth() / treeLineages.getMaxWidth();
            scaleY = getBounds().getHeight() / (treeLineages.getMaxHeight() - 1.0);
//            isCalibrated = true;
        }

        ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        float count = 0.0F;
        for (TreeLineages.Lineage root : treeLineages.getRootLineages()) {
            float hue = count / treeLineages.getRootLineages().size();
            paintTreeLineages((Graphics2D)g, hue, root, treeLineages.getMaxWidth(), 0.0);

            count += 1.0F;
        }
    }

    private void paintTreeLineages(Graphics2D g2, float hue, TreeLineages.Lineage lineage, double x, double y) {
        double x1 = x + lineage.dx;
        double y1 = y + lineage.dy;

        paintLineage(g2, hue, lineage.child1, x1, y1);
        paintLineage(g2, hue, lineage.child2, x1, y1);
    }

    private void paintLineage(Graphics2D g2, float hue, TreeLineages.Lineage lineage, double x, double y) {
        double x1 = x + lineage.dx;
        double y1 = y + lineage.dy;

//        g2.setColor(Color.getHSBColor(hue, 1.0F, 1.0F));
        g2.setColor(new Color(0.0F, 0.0F, 0.0F, 0.1F));
        g2.setStroke(new BasicStroke(0.5F));
        g2.draw(new Line2D.Double(transformX(x), transformY(y), transformX(x1), transformY(y1)));

        if (lineage.child1 != null) {
            paintLineage(g2, hue, lineage.child1, x1, y1);
            paintLineage(g2, hue, lineage.child2, x1, y1);
        } else {
            paintTip(g2, lineage, x1, y1);
        }
    }

    private void paintTip(Graphics2D g2, TreeLineages.Lineage lineage, double x, double y) {
        Shape node = new Ellipse2D.Float((float)transformX(x) - nodeWidth * 0.5F, (float)transformY(y) - nodeWidth * 0.5F, nodeWidth, nodeWidth);
        g2.setColor(Color.getHSBColor(((float)lineage.tipNumber) / (float)treeLineages.getMaxHeight(), 0.5F, 1.0F));
        g2.fill(node);
        g2.setColor(Color.black);
        g2.draw(node);
    }

    private double transformX(double x) {
        return (x * scaleX);
    }

    private double transformY(double y) {
        return (y * scaleY);
    }

    private TreeLineages treeLineages;

    private float xLabelOffset = 5.0F;
    private float yLabelOffset = 5.0F;
    private boolean isCalibrated = false;

    private double scaleX = 1.0;
    private double scaleY = 1.0;

    private float nodeWidth = 10.0F;
}
