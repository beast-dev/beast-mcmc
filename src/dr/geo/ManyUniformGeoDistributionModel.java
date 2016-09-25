/*
 * ManyUniformGeoDistributionModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.geo;

import dr.inference.model.*;

import java.util.*;

/**
 * @author Marc A. Suchard
 */
public class ManyUniformGeoDistributionModel extends AbstractModelLikelihood {


    public ManyUniformGeoDistributionModel(String name,
                                           List<Parameter> parameters,
                                           List<GeoSpatialDistribution> distributions) {
        this(name, parameters, distributions, null);
    }

    public ManyUniformGeoDistributionModel(String name,
                                           List<Parameter> parameters,
                                           List<GeoSpatialDistribution> distributions,
                                           List<Likelihood> oldLikelihoods) {
        super(name);

        assert (parameters.size() == distributions.size());
        count = parameters.size();

        this.parameters = parameters;
        this.distributions = distributions;
        this.oldLikelihoods = oldLikelihoods;

        parameterMap = new HashMap<Parameter, Integer>(count);
        likelihood = new double[count];
        storedLikelihood = new double[count];
        likelihoodKnown = new boolean[count];
        storedLikelihoodKnown  = new boolean[count];

        for (int i = 0; i < count; ++i) {
            final Parameter parameter = parameters.get(i);
            addVariable(parameter);
            parameterMap.put(parameter, i);
        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public final double getLogLikelihood() {
        if (allLikelihoodsKnown) {
            logLikelihood = calculateLogLikelihood();
            allLikelihoodsKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        allLikelihoodsKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        if (oldLikelihoods == null) {
            return new HashSet<Likelihood>(Arrays.asList(this));
        } else {
            Set<Likelihood> likelihoodSet = new HashSet<Likelihood>(oldLikelihoods);
            likelihoodSet.add(this);
            return likelihoodSet;
        }
    }

    @Override
    protected void storeState() {
        storedAllLikelihoodsKnown = allLikelihoodsKnown;
        storedLogLikelihood = logLikelihood;

        System.arraycopy(likelihood, 0, storedLikelihood, 0, count);
        System.arraycopy(likelihoodKnown, 0, storedLikelihoodKnown, 0 , count);
    }

    @Override
    protected void restoreState() {
        allLikelihoodsKnown = storedAllLikelihoodsKnown;
        logLikelihood = storedLogLikelihood;

        double[] tmp1 = likelihood;
        likelihood = storedLikelihood;
        storedLikelihood = tmp1;

        boolean[] tmp2 = likelihoodKnown;
        likelihoodKnown = storedLikelihoodKnown;
        storedLikelihoodKnown = tmp2;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        assert (variable instanceof Parameter);

        int whichParameter = parameterMap.get(variable);

        likelihoodKnown[whichParameter] = false;
        allLikelihoodsKnown = false;
    }

    private double calculateLogLikelihood() {
        double logLikelihood = 0;

        for (int i = 0; i < count; ++i) {
            if (!likelihoodKnown[i]) {
                likelihood[i] = distributions.get(i).logPdf(
                        parameters.get(i).getParameterValues()
                );
            }
            logLikelihood += likelihood[i];
        }

        return logLikelihood;
    }

    private final List<Parameter> parameters;
    private final List<GeoSpatialDistribution> distributions;
    private final List<Likelihood> oldLikelihoods;

    private final int count;

    private final Map<Parameter, Integer> parameterMap;

    private boolean[] likelihoodKnown;
    private boolean[] storedLikelihoodKnown;

    private double[] likelihood;
    private double[] storedLikelihood;

    private boolean allLikelihoodsKnown = false;
    private boolean storedAllLikelihoodsKnown;

    private double logLikelihood;
    private double storedLogLikelihood;
}
