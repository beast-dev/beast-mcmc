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
            scaleX = ((double)getWidth()) / treeLineages.getMaxWidth();
            scaleY = ((double)getHeight()) / (treeLineages.getMaxHeight() + 10.0);

//            isCalibrated = true;
        }

        ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        for (TreeLineages.Lineage root : treeLineages.getRootLineages()) {
            paintLineage((Graphics2D)g, root, treeLineages.getMaxWidth(), 0.0);
        }
    }

    private void paintLineage(Graphics2D g2, TreeLineages.Lineage lineage, double x, double y) {
        double x1 = x + lineage.dx;
        double y1 = y + lineage.dy;

        if (lineage.child1 != null) {
            paintLineage(g2, lineage.child1, x1, y1);
            paintLineage(g2, lineage.child2, x1, y1);
        }

        g2.setColor(Color.black);
        g2.setStroke(new BasicStroke(0.25F));
        g2.draw(new Line2D.Double(transformX(x), transformY(y), transformX(x1), transformY(y1)));
    }

//    private void paintTip(Graphics2D g2, Clade clade) {
//        float x = transformX(clade.x);
//        float y = transformY(clade.y);
//
//        Shape node = new Ellipse2D.Float(x, y, (float)(nodeWidth * scaleX), (float)(nodeHeight * scaleY));
//        g2.setColor(Color.yellow);
//        g2.fill(node);
//        g2.setColor(Color.black);
//        g2.draw(node);
//        if (clade.label != null) {
//            g2.drawString(clade.label, x + xLabelOffset, yLabelOffset);
//        }
//    }

    private double transformX(double x) {
        return (x * scaleX);
    }

    private double transformY(double y) {
        return (y * scaleY);
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
