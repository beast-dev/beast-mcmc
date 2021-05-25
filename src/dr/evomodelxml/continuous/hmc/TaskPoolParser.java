/*
 * FullyConjugateTreeTipsPotentialDerivativeParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.util.TaskPool;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

@SuppressWarnings("unused")
public class TaskPoolParser extends AbstractXMLObjectParser {

    private static final String TAXON_PARSER_NAME = "taxonTaskPool";
    private static final String TASk_PARSER_NAME = "taskPool";
    public static final String THREAD_COUNT = "threadCount";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        GradientWrtParameterProvider gradient = (GradientWrtParameterProvider)
                xo.getChild(GradientWrtParameterProvider.class);

        int taskCount = (tree != null) ? tree.getExternalNodeCount() : gradient.getDimension();
        int threadCount = xo.getAttribute(THREAD_COUNT, 1);

        return new TaskPool(taskCount, threadCount);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "A thread pool for per-taxon specific operations";
    }

    @Override
    public Class getReturnType() {
        return TaskPool.class;
    }

    @Override
    public String getParserName() {
        return TAXON_PARSER_NAME;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { TAXON_PARSER_NAME, TASk_PARSER_NAME };
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(THREAD_COUNT, true),
            new XORRule(
                    new ElementRule(Tree.class),
                    new ElementRule(GradientWrtParameterProvider.class)
            ),
    };
}
