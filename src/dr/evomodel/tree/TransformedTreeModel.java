/*
 * TransformedTreeModel.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * A transformation of a tree model
 *
 * @author Marc Suchard
 */
public class TransformedTreeModel extends TreeModel implements Citable {

    public TransformedTreeModel(String id, Tree tree, Parameter scale) {
        super(id, tree);

        scale.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        this.scale = scale;
        addVariable(scale);

        Logger log = Logger.getLogger("dr.evomodel.tree");
        log.info("Creating a transform tree. Linear scaling parameter: " + scale.getId() + "\n\tPlease cite:");
        log.info(Citable.Utils.getCitationString(this));
    }

    public double getNodeHeight(NodeRef node) {
        return transform(node, super.getNodeHeight(node));
    }

    private double transform(NodeRef node, double originalHeight) {
        if (isExternal(node)) {
            return originalHeight;
        }
        final double rootHeight = getRootHeightParameter().getParameterValue(0);
        return rootHeight - scale.getParameterValue(0) * (rootHeight - originalHeight);
    }


    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == scale) {
            pushTreeChangedEvent(); // All internal node heights have changed!
        } else {
            super.handleVariableChangedEvent(variable, index, type);
        }
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                CommonCitations.LEMEY_MIXTURE_2012
        );
        return citations;
    }

    Parameter scale;
}
