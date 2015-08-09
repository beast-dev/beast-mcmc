/*
 * AbstractDistributionLikelihood.java
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

package dr.inference.distribution;

import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.util.Attribute;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: AbstractDistributionLikelihood.java,v 1.4 2005/05/24 20:25:59 rambaut Exp $
 */

public abstract class AbstractDistributionLikelihood extends Likelihood.Abstract {

    public AbstractDistributionLikelihood(Model model) {

        super(model);
    }

    /**
     * Adds a statistic, this is the data for which the likelihood is calculated.
     *
     * @param data to add
     */
    public void addData(Attribute<double[]> data) {
        dataList.add(data);
    }

    protected ArrayList<Attribute<double[]>> dataList = new ArrayList<Attribute<double[]>>();

    public abstract double calculateLogLikelihood();

    /**
     * Overridden to always return false.
     */
    protected boolean getLikelihoodKnown() {
        return false;
	}

    public List<Attribute<double[]>> getDataList() {
        return dataList;
    }
}

