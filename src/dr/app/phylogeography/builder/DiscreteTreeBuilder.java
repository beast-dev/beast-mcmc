package dr.app.phylogeography.builder;

import dr.app.phylogeography.structure.Coordinates;
import dr.app.phylogeography.structure.Layer;
import dr.app.phylogeography.structure.Line;
import dr.app.phylogeography.structure.Style;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import jam.panels.OptionsPanel;
import org.virion.jam.components.RealNumberField;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTreeBuilder extends AbstractBuilder {

    private static final String BUILDER_NAME = "Discrete Diffusion Tree";

    public static final String LONGITUDE_ATTRIBUTE = "long";
    public static final String LATITUDE_ATTRIBUTE = "lat";
    public static final String TIME_ATTRIBUTE = "height";

    private JPanel editPanel = null;

    private final RealNumberField maxAltitudeField = new RealNumberField(0, Double.MAX_VALUE);

    private double maxAltitude = 0;

    public DiscreteTreeBuilder() {
    }

    protected Layer buildLayer() throws BuildException {
        Layer layer = new Layer(getName(), getDescription(), isVisible());
        buildTree(layer, getDataFile().getFirstTree());

        return layer;
    }

    public JPanel getEditPanel() {
        if (editPanel == null) {
            OptionsPanel editPanel = new OptionsPanel();
            maxAltitudeField.setColumns(10);
            editPanel.addComponentWithLabel("Maximum Altitude:", maxAltitudeField);
            this.editPanel = editPanel;
        }
        maxAltitudeField.setValue(maxAltitude);
        return editPanel;
    }

    public void setFromEditPanel() {
        maxAltitude = maxAltitudeField.getValue(0.0);
        invalidate();
    }

    private void buildTree(Layer layer, final Tree tree) throws BuildException {
        buildTree(layer, tree, tree.getRoot());
    }

    private void buildTree(Layer layer, final Tree tree, final NodeRef node) throws BuildException {
        if (!tree.isRoot(node)) {
            NodeRef parent = tree.getParent(node);
            double long0 = getDoubleAttribute(tree, parent, LONGITUDE_ATTRIBUTE);
            double lat0 = getDoubleAttribute(tree, parent, LATITUDE_ATTRIBUTE);
            double time0 = getDoubleAttribute(tree, parent, TIME_ATTRIBUTE);
            Style style0 = new Style(Color.red, 1.0);

            double long1 = getDoubleAttribute(tree, node, LONGITUDE_ATTRIBUTE);
            double lat1 = getDoubleAttribute(tree, node, LATITUDE_ATTRIBUTE);
            double time1 = getDoubleAttribute(tree, node, TIME_ATTRIBUTE);
            Style style1 = new Style(Color.red, 1.0);

            double duration = -1.0;

            Line line = new Line(
                    "",
                    new Coordinates(long0, lat0), time0, style0,
                    new Coordinates(long1, lat1), time1, style1,
                    maxAltitude,
                    duration
            );
            layer.addItem(line);
        }
        if (!tree.isExternal(node)) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                buildTree(layer, tree, child);
            }
        }
    }

    private double getDoubleAttribute(Tree tree, NodeRef node, String attributeName) throws BuildException {
        Double value = (Double)tree.getNodeAttribute(node, attributeName);
        if (value == null) {
            throw new BuildException("Tree doesn't have attribute, " + attributeName + ", for one or more nodes");
        }
        return value;
    }
    public String getBuilderName() {
        return FACTORY.getBuilderName();
    }

    public final static BuilderFactory FACTORY = new BuilderFactory() {

        public String getBuilderName() {
            return BUILDER_NAME;
        }

        public Builder createBuilder() {
            return new DiscreteTreeBuilder();
        }
    };

}
