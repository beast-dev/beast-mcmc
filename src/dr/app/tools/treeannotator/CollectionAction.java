/*
 * CollectionAction.java
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

package dr.app.tools.treeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

class CollectionAction implements CladeAction {
    private final Set<String> attributeNames = new LinkedHashSet<>();

    public void addAttributeName(String attributeName) {
        this.attributeNames.add(attributeName);
    }
    public void addAttributeNames(Collection<String> attributeNames) {
        this.attributeNames.addAll(attributeNames);
    }

    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        collectAttributesForClade(clade, attributeNames, tree, node);
    }

    @Override
    public boolean expectAllClades() {
        return true;
    }

    private void collectAttributesForClade(Clade clade, Set<String> attributeNames, Tree tree, NodeRef node) {
        int i = 0;
        Object[] values = new Object[attributeNames.size()];
        for (String attributeName : attributeNames) {
            Object value;
            if (attributeName.equals("height")) {
                value = tree.getNodeHeight(node);
                ((BiClade)clade).addHeightValue((Double)value);
            } else if (attributeName.equals("length")) {
                value = tree.getBranchLength(node);
            } else {
                value = tree.getNodeAttribute(node, attributeName);
                if (value instanceof String && ((String) value).startsWith("\"")) {
                    value = ((String) value).replaceAll("\"", "");
                }
            }

            values[i] = value;
            i++;
        }
        clade.addAttributeValues(values);
    }
}
