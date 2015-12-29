/*
 * MixedDistributionLikelihood.java
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

import dr.inference.model.CompoundModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Statistic;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that returns the log likelihood of a multidimensional statistic
 * being distributed according to a mixture of parametric distributions, where
 * the distribution membership of each element is determined by a separate vector
 * of indices.
 *
 * @author Alexei Drummond
 * @version $Id: DistributionLikelihood.java,v 1.11 2005/05/25 09:35:28 rambaut Exp $
 */

public class MixedDistributionLikelihood extends Likelihood.Abstract {

    public MixedDistributionLikelihood(ParametricDistributionModel[] distributions, Statistic data, Statistic indicators) {
        // Mixed Distribution Likelihood contains two distribution models, not necessarily constant.
        // To cater for that, they need to be returned as the "Model" of this likelyhood so that their state is correctly
        // restored when an operation involving their parameters fails.

        super(new CompoundModel("MixedDistributions"));

        final CompoundModel cm = (CompoundModel) this.getModel();
        for (ParametricDistributionModel m : distributions) {
            cm.addModel(m);
        }

        this.distributions = distributions;
        this.data = data;
        this.indicators = indicators;

        if (indicators.getDimension() == data.getDimension() - 1) {
            impliedOne = true;
        } else if (indicators.getDimension() != data.getDimension()) {
            throw new IllegalArgumentException("Indicators length (" + indicators.getDimension() +
                    ") != data length (" + data.getDimension() + ")");
        }
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public final double calculateLogLikelihood() {

        double logL = 0.0;

        for (int j = 0; j < data.getDimension(); j++) {

            int index;
            if (impliedOne) {
                if (j == 0) {
                    index = 1;
                } else {
                    index = (int) indicators.getStatisticValue(j - 1);
                }
            } else {
                index = (int) indicators.getStatisticValue(j);
            }

            logL += distributions[index].logPdf(data.getStatisticValue(j));
        }
        //System.err.println("mixed: " + logL);
        return logL;
    }

    /**
     * Overridden to always return false.
     */
    protected boolean getLikelihoodKnown() {
        return false;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    private final ParametricDistributionModel[] distributions;
    private final Statistic data;
    private final Statistic indicators;
    private boolean impliedOne = false;

    public Model[] getUniqueModels() {
        Model[] m = new Model[distributions[0] == distributions[1] ? 1 : 2];
        m[0] = distributions[0];
        if (m.length > 1) m[1] = distributions[1];
        return m;
    }
}
