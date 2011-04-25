package dr.app.treespace;

import dr.app.treespace.CladeSystem.Clade;
import jebl.evolution.graphs.Node;
import jebl.evolution.trees.RootedTree;
import jebl.math.Random;
import sun.awt.geom.Curve;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.image.ColorModel;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreePlotter extends JComponent {

    private final static Color[] COLOUR_SCHEME = {
            new Color(0.879F, 0.261F, 0.262F),	//Africa,
            new Color(0.816F, 0.765F, 0.376F),	//USA,
            new Color(0.26F, 0.141F, 0.632F),	//Taiwan,
            new Color(0.242F, 0.445F, 0.718F),	//China,
            new Color(0.592F, 0.669F, 0.295F),	//Russia,
            new Color(0.791F, 0.27F, 0.146F),	//Oceania,
            new Color(0.359F, 0.425F, 0.833F),	//Asia,
            new Color(0.374F, 0.623F, 0.505F),	//Japan,
            new Color(0.785F, 0.585F, 0.209F),	//Mexico,
            new Color(0.917F, 0.58F, 0.322F),	//South America,
            new Color(0.64F, 0.46F, 0.28F),	    //Canada,
            new Color(0.599F, 0.772F, 0.513F),	//Europe,
            new Color(0.551F, 0.242F, 0.598F),	//Southeast Asia,
            new Color(0.43F, 0.674F, 0.744F)	//South Korea
    };

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
            scaleX = (getBounds().getWidth() - borderWidth - borderWidth) / treeLineages.getMaxWidth();
            scaleY = (getBounds().getHeight() - borderWidth - borderWidth) / (treeLineages.getMaxHeight() - 1.0);
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
        double y1 = lineage.dy;

        paintLineage(g2, hue, lineage.child1, x1, y1);
        paintLineage(g2, hue, lineage.child2, x1, y1);
    }

    private void paintLineage(Graphics2D g2, float hue, TreeLineages.Lineage lineage, double x, double y) {
        double x1 = x + lineage.dx;
        double y1 = lineage.dy;

        if (lineage.child1 != null) {
            paintLineage(g2, hue, lineage.child1, x1, y1);
            paintLineage(g2, hue, lineage.child2, x1, y1);
        } else {
            paintTip(g2, lineage, x1, y1);
        }

//        g2.setColor(Color.getHSBColor(hue, 1.0F, 1.0F));
//        g2.setColor(new Color(0.0F, 0.0F, 0.0F, 0.1F));
//        g2.setColor(COLOUR_SCHEME[lineage.state]);
        g2.setColor(new Color(COLOUR_SCHEME[lineage.state].getRed(), COLOUR_SCHEME[lineage.state].getGreen(), COLOUR_SCHEME[lineage.state].getBlue(), 32));

        float thickness = 0.25F;

        thickness = 2.5F * (Math.min((float)lineage.tipCount, 10.0F) / 10.0F);

        g2.setStroke(new BasicStroke(thickness));
//        Shape curve = new Line2D.Double(transformX(x), transformY(y), transformX(x1), transformY(y1));

        double xp1 = transformX(x);
        double yp1 = transformY(y);
        double xp2 = transformX(x1);
        double yp2 = transformY(y1);
        double xc = (xp1 + xp2) / 2.0;

        Shape curve = new CubicCurve2D.Double(xp1, yp1, xc, yp1, xc, yp2, xp2, yp2);

        g2.draw(curve);

    }

    private void paintTip(Graphics2D g2, TreeLineages.Lineage lineage, double x, double y) {
        Shape node = new Ellipse2D.Float((float)transformX(x) - nodeWidth * 0.5F, (float)transformY(y) - nodeWidth * 0.5F, nodeWidth, nodeWidth);
//        g2.setColor(Color.getHSBColor(((float)lineage.tipNumber) / (float)treeLineages.getMaxHeight(), 0.5F, 1.0F));
//        g2.setColor(COLOUR_SCHEME[lineage.state]);
        g2.setColor(new Color(COLOUR_SCHEME[lineage.state].getRed(), COLOUR_SCHEME[lineage.state].getGreen(), COLOUR_SCHEME[lineage.state].getBlue(), 32));
        g2.fill(node);
//        g2.setColor(Color.black);
//        g2.draw(node);
    }

    private double transformX(double x) {
        return (x * scaleX) + borderWidth;
    }

    private double transformY(double y) {
        return (y * scaleY) + borderWidth;
    }

    private TreeLineages treeLineages;

    private float xLabelOffset = 5.0F;
    private float yLabelOffset = 5.0F;
    private boolean isCalibrated = false;

    private double scaleX = 1.0;
    private double scaleY = 1.0;

    private float nodeWidth = 5.0F;
    private float borderWidth = 10.0F;
}
