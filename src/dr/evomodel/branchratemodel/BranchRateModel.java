/*
 * BranchRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.*;
import dr.inference.model.Model;

/**
 * Date: Dec 13, 2004
 * Time: 1:59:24 PM
 *
 * @author Alexei Drummond
 * @version $Id: BranchRateModel.java,v 1.4 2005/05/24 20:25:57 rambaut Exp $
 */
public interface BranchRateModel extends Model, BranchRates, TreeTraitProvider, TreeTrait<Double> {
    public static final String BRANCH_RATES = "branchRates";
    public static final String RATE = "rate";

    // This is inherited from BranchRates:
    // double getBranchRate(Tree tree, NodeRef node);
}
