package dr.app.phylogeography.builder;

import dr.app.phylogeography.structure.*;
import dr.app.phylogeography.spread.SpreadDocument;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.graphs.Node;

import javax.swing.*;
import java.awt.*;

import jam.panels.OptionsPanel;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.components.RealNumberField;

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

    public Layer buildLayer() {
        Layer layer = new Layer(getName(), getDescription(), isVisible());
        buildTree(layer, (RootedTree)getDataFile().getFirstTree());
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
    }

    private void buildTree(Layer layer, final RootedTree tree) {
        buildTree(layer, tree, tree.getRootNode());
    }

    private void buildTree(Layer layer, final RootedTree tree, final Node node) {
        if (!tree.isRoot(node)) {
            Node parent = tree.getParent(node);
            double long0 = (Double)parent.getAttribute(LONGITUDE_ATTRIBUTE);
            double lat0 = (Double)parent.getAttribute(LATITUDE_ATTRIBUTE);
            double time0 = (Double)parent.getAttribute(TIME_ATTRIBUTE);
            Style style0 = new Style(Color.red, 1.0);

            double long1 = (Double)node.getAttribute(LONGITUDE_ATTRIBUTE);
            double lat1 = (Double)node.getAttribute(LATITUDE_ATTRIBUTE);
            double time1 = (Double)node.getAttribute(TIME_ATTRIBUTE);
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
        for (Node child : tree.getChildren(node)) {
            buildTree(layer, tree, child);
        }
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
