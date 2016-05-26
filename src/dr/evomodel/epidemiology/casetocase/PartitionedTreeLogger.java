/*
 * PartitionedTreeLogger.java
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

package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.inference.loggers.LogFormatter;

import java.text.NumberFormat;

/**
 * Logs a partitioned tree. Breaks branches in two at infection times and adds a node of degree 2.
 * TreeAnnotator won't work.
 *
 * @author Matthew Hall
 */
public class PartitionedTreeLogger extends TreeLogger {

    CaseToCaseTreeLikelihood c2cTL;
    Tree originalTree;
    private LogUpon condition = null;
    private TreeTraitProvider.Helper treeTraitHelper;

    public PartitionedTreeLogger(CaseToCaseTreeLikelihood c2cTL, Tree tree, BranchRates branchRates,
                      TreeAttributeProvider[] treeAttributeProviders,
                      TreeTraitProvider[] treeTraitProviders,
                      LogFormatter formatter, int logEvery, boolean nexusFormat,
                      boolean sortTranslationTable, boolean mapNames, NumberFormat format,
                      TreeLogger.LogUpon condition) {
        super(tree, branchRates, treeAttributeProviders, treeTraitProviders, formatter, logEvery, nexusFormat,
                sortTranslationTable, mapNames, format, condition);
        this.c2cTL = c2cTL;
        this.originalTree = tree;
        this.condition = condition;
    }

    public void log(long state){

        final boolean doIt = condition != null ? condition.logNow(state) :
                (logEvery < 0 || ((state % logEvery) == 0));
        if(!doIt)
            return;

        setTree(c2cTL.addTransmissionNodes(originalTree));

        super.log(state);

    }



}
