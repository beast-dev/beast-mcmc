package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockType;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeModelParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;

/**
 * @author Alexei Drummond
 */
public class TreeModelGenerator extends Generator {

    public TreeModelGenerator(BeautiOptions options) {
        super(options);
    }

    /**
     * Write tree model XML block.
     *
     * @param writer the writer
     */
    void writeTreeModel(XMLWriter writer) {

        writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>("id", "treeModel"), false);

        if (options.userTree) {
            writer.writeTag("tree", new Attribute.Default<String>("idref", "startingTree"), true);
        } else {
            writer.writeTag(CoalescentSimulator.COALESCENT_TREE, new Attribute.Default<String>("idref", "startingTree"), true);
        }

        writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("id", "treeModel.rootHeight"), true);
        writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);


        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS, new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"));
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("id", "treeModel.internalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
                new Attribute[]{
                        new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                        new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true")
                });
        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("id", "treeModel.allInternalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        if (options.clockType == ClockType.RANDOM_LOCAL_CLOCK) {
            writer.writeOpenTag(TreeModelParser.NODE_RATES,
                    new Attribute[]{
                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("id", "localClock.rates"), true);
            writer.writeCloseTag(TreeModelParser.NODE_RATES);

            writer.writeOpenTag(TreeModelParser.NODE_TRAITS,
                    new Attribute[]{
                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>("id", "localClock.changes"), true);
            writer.writeCloseTag(TreeModelParser.NODE_TRAITS);
        }

        writer.writeCloseTag(TreeModel.TREE_MODEL);
    }
}
