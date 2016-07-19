/*
 * SingleTipObservationProcess.java
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

package dr.oldevomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteRateModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodelxml.MSSD.SingleTipObservationProcessParser;
import dr.inference.model.Parameter;

/**
 * Package: SingleTipObservationProcess
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 19, 2008
 * Time: 2:57:14 PM
 */
public class SingleTipObservationProcess extends AnyTipObservationProcess {
    protected Taxon sourceTaxon;

    public SingleTipObservationProcess(TreeModel treeModel, PatternList patterns, SiteRateModel siteModel,
                                       BranchRateModel branchRateModel, Parameter mu, Parameter lam, Taxon sourceTaxon) {
        super(SingleTipObservationProcessParser.MODEL_NAME, treeModel, patterns, siteModel, branchRateModel, mu, lam);
        this.sourceTaxon = sourceTaxon;
    }

    public double calculateLogTreeWeight() {
        return -lam.getParameterValue(0) / (getAverageRate() * mu.getParameterValue(0));
    }

}
