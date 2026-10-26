/*
 * SampledTipHeights.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.NewTreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;

import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Finds sampled tip heights (tip-date sampling: a leaf height handed out as a parameter) among the objects parsed so
 * far, so that parsers can refuse analytic gradients with respect to them: the tree likelihood and tree prior terms in
 * the tip heights are not differentiated yet. The trees are looked up in the parser's store, so the tree models need
 * no bookkeeping of their own.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SampledTipHeights {

    public static final String GRADIENT_NOT_IMPLEMENTED =
            "Gradients with respect to sampled tip dates are not yet implemented";

    private SampledTipHeights() {
    }

    /**
     * @param store the objects parsed so far (AbstractXMLObjectParser.getStore())
     * @return true if the parameter is, or contains (as a compound parameter), the height of a leaf of a parsed tree
     * that samples its tip dates
     */
    public static boolean containsSampledTipHeight(Parameter parameter, Map<String, XMLObject> store) {
        if (parameter == null || store == null) {
            return false;
        }
        Set<Object> trees = Collections.newSetFromMap(new IdentityHashMap<Object, Boolean>());
        for (XMLObject xo : store.values()) {
            Object object = xo.getNativeObject();
            if ((object instanceof DefaultTreeModel && ((DefaultTreeModel) object).isTipDateSampled()) ||
                    (object instanceof NewTreeModel && ((NewTreeModel) object).isTipDateSampled())) {
                trees.add(object);
            }
        }
        return !trees.isEmpty() && contains(parameter, trees);
    }

    private static boolean contains(Parameter parameter, Set<Object> trees) {
        if (parameter instanceof CompoundParameter) {
            CompoundParameter compound = (CompoundParameter) parameter;
            for (int i = 0; i < compound.getParameterCount(); i++) {
                if (contains(compound.getParameter(i), trees)) {
                    return true;
                }
            }
            return false;
        }
        for (Object tree : trees) {
            if (tree instanceof DefaultTreeModel) {
                DefaultTreeModel treeModel = (DefaultTreeModel) tree;
                NodeRef node = treeModel.getNodeOfParameter(parameter);
                if (node != null && treeModel.isExternal(node)) {
                    return true;
                }
            } else {
                NewTreeModel treeModel = (NewTreeModel) tree;
                NodeRef node = treeModel.getNodeOfParameter(parameter);
                if (node != null && treeModel.isExternal(node)) {
                    return true;
                }
            }
        }
        return false;
    }
}
