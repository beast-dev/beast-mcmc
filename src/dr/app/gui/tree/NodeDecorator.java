/*
 * NodeDecorator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.gui.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

import java.awt.*;

/**
 * @author Alexei Drummond
 * Date: Dec 3, 2004
 * Time: 9:57:30 AM
 *
 * @version $Id: NodeDecorator.java,v 1.1 2005/11/29 12:59:36 rambaut Exp $
 */
public interface NodeDecorator {

    boolean isDecoratable(Tree tree, NodeRef node);

    void decorateNode(Tree tree, NodeRef node, Graphics2D g2, CoordinateTransform transform);

}
