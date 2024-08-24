/*
 * TreePartitionData.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
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
    public TreePartitionData(BeautiOptions options, String name, String fileName, TreeHolder trees) {
        super(options, name, fileName);

        this.trees = trees;
    }

    @Override
    public String getPrefix() {
        return getName() + ".";
    }

    @Override
    public TaxonList getTaxonList() {
        return trees.getTrees().get(0);
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

    public TreeHolder getTrees() {
        return trees;
    }

    public int getTreeCount() {
        return trees.getTreeCount();
    }


    private final TreeHolder trees;
}
