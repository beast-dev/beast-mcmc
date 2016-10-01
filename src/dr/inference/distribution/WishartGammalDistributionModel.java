/*
 * WishartGammalDistributionModel.java
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

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.WishartDistribution;
import dr.math.distributions.WishartStatistics;

/**
 * A class that acts as a model for multivariate Wishart-scaled-by-gamma distribution
 * See:
 *  Huang A and Wand MP (in press) Simple marginally noninformative prior distributions for covariance matrices.
 *  Bayesian Analysis.
 *
 * @author Marc Suchard
 * @author Max Tolkoff
 * @author Philipe Lemey
 * @author Gabriela Cybis
 */

public class WishartGammalDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel,
        WishartStatistics {

    public WishartGammalDistributionModel(Parameter df, Parameter mixing, Parameter scale, boolean randomMixing) {
        super(WishartGammalDistributionModel.TYPE);

        this.df = df;
        this.mixing = mixing;
        this.scale = scale;

        addVariable(df);
        addVariable(mixing);
        addVariable(scale);

        this.randomMixing = randomMixing;

        df.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                df.getDimension()));

        mixing.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                mixing.getDimension()));

        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                scale.getDimension()));

        dim = mixing.getDimension();
        if (dim != scale.getDimension()) {
            throw new IllegalArgumentException();
        }

        wishartDistributionKnown = false;
    }

    // *****************************************************************
    // Interface MultivariateDistribution
    // *****************************************************************


    public double logPdf(double[] x) {
        if (x.length != dim * dim) {
            throw new IllegalArgumentException("Data of wrong dimension in " + getId());
        }

        if (!wishartDistributionKnown) {
            wishart = createNewWishartDistribution();
            wishartDistributionKnown = true;
        }

        double logLike = wishart.logPdf(x);
        if (randomMixing) {
            for (int i = 0; i < dim; ++i) {
                double mixing = this.mixing.getParameterValue(i);
                double scale = this.scale.getParameterValue(i);
                logLike += GammaDistribution.logPdf(mixing, 0.5, scale * scale);
            }
        }
        return logLike;
    }

    public double[][] getScaleMatrix() {
        double df = this.df.getParameterValue(0);
        double[][] scaleMatrix = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            scaleMatrix[i][i] = 1.0 / (2.0 * df * mixing.getParameterValue(i));
        }
        return scaleMatrix;
    }

    public double getDF() {
        return df.getParameterValue(0) + dim - 1;
    }

    public double[] getMean() {
        return new double[dim];
    }

    public String getType() {
        return TYPE;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == mixing || variable == df) {
            wishartDistributionKnown = false;
        }

        if (!randomMixing && variable == mixing) {
            throw new IllegalArgumentException("WishartGamma distribution is not configured for random mixing");
        }
    }

    protected void storeState() {
        storedWishart = wishart;
        storedWishartDistributionKnown = wishartDistributionKnown;
    }

    protected void restoreState() {
        wishartDistributionKnown = storedWishartDistributionKnown;
        wishart = storedWishart;
    }

    protected void acceptState() {
    } // no additional state needs accepting


    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        throw new UnsupportedOperationException("Not implemented");
    }

    // **************************************************************
    // Private instance variables and functions
    // **************************************************************

    private WishartDistribution createNewWishartDistribution() {
        return new WishartDistribution(getDF(), getScaleMatrix());
    }

    private final int dim;

    private final Parameter df;
    private final Parameter mixing;
    private final Parameter scale;

    private final boolean randomMixing;

    private WishartDistribution wishart;
    private WishartDistribution storedWishart;

    private boolean wishartDistributionKnown;
    private boolean storedWishartDistributionKnown;

    public static final String TYPE = "WishartGamma";

	@Override
	public double[] nextRandom() {
		// TODO Auto-generated method stub
		return null;
	}
}
