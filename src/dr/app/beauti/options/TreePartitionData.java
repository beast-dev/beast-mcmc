/*
 * Copyright (c) 2024. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
 * Morbi non lorem porttitor neque feugiat blandit. Ut vitae ipsum eget quam lacinia accumsan.
 * Etiam sed turpis ac ipsum condimentum fringilla. Maecenas magna.
 * Proin dapibus sapien vel ante. Aliquam erat volutpat. Pellentesque sagittis ligula eget metus.
 * Vestibulum commodo. Ut rhoncus gravida arcu.
 */

package dr.app.beauti.options;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TreeDataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;

/**
 * A concrete implementation of AbstractPartitionData that provides a tree-as-data partition for
 * use with models like 'ThorneyBEAST'
 * @author Andrew Rambaut
 */
public class TreePartitionData extends AbstractPartitionData {
    public TreePartitionData(BeautiOptions options, String name, String fileName, Tree tree) {
        super(options, name, fileName);

        this.tree = tree;
    }

    @Override
    public String getPrefix() {
        return getName() + ".";
    }

    @Override
    public TaxonList getTaxonList() {
        return tree;
    }

    @Override
    public int getSiteCount() {
        return 0;
    }

    @Override
    public int getPatternCount() {
        return 0;
    }

    @Override
    public DataType getDataType() {
        return TreeDataType.INSTANCE;
    }

    @Override
    public String getDataDescription() {
        return TreeDataType.INSTANCE.getDescription();
    }

    private final Tree tree;
}
