/*
 * RLTVLogger.java
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

package dr.evomodel.tree.randomlocalmodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

/**
 * A logger for a random local tree variable.
 * It logs only the selected parameters in pairs of node number, variable value
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class RLTVLogger extends MCLogger {

    private TreeModel treeModel;
    private RandomLocalTreeVariable randomLocal;

    public RLTVLogger(LogFormatter formatter,
                      int logEvery, TreeModel treeModel,
                      RandomLocalTreeVariable randomLocal) {

        super(formatter, logEvery, false);
        this.treeModel = treeModel;
        this.randomLocal = randomLocal;
    }

    public void startLogging() {
        logLine("State\tRate changes");
    }

    public void log(long state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {

            int nodeCount = treeModel.getNodeCount();

            StringBuilder builder = new StringBuilder();
            builder.append(state);
            for (int i = 0; i < nodeCount; i++) {

                NodeRef node = treeModel.getNode(i);

                if (randomLocal.isVariableSelected(treeModel, node)) {
                    builder.append("\t");
                    builder.append(node.getNumber());
                    builder.append("\t");
                    builder.append(randomLocal.getVariable(treeModel, node));
                }
            }

            logLine(builder.toString());
        }
    }

    public void stopLogging() {

        super.stopLogging();
    }

}
