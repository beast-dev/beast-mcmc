/*
 * ColourSampler.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evolution.colouring;

import dr.evolution.coalescent.structure.MetaPopulation;
import dr.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 * @author Gerton Lunter
 * @author Andrew Rambaut
 * @version $Id: ColourSampler.java,v 1.5 2006/09/11 09:33:01 gerton Exp $
 */
public interface ColourSampler {

    DefaultTreeColouring sampleTreeColouring(Tree tree, ColourChangeMatrix colourChangeMatrix, MetaPopulation mp);

    double getProposalProbability(TreeColouring treeColouring, Tree tree, ColourChangeMatrix colourChangeMatrix, MetaPopulation mp);

    int[] getLeafColourCounts();
}