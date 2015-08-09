/*
 * SpeciesTreeSimplePrior.java
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

package dr.evomodel.speciation;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.NodeRef;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.CompoundModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Joseph Heled
 *         Date: 8/12/2008
 */
public class SpeciesTreeSimplePrior extends Likelihood.Abstract {

    private final SpeciesTreeModel sTree;
    //private final ParametricDistributionModel dist;
    private final ParametricDistributionModel tips;
    private final Parameter sigma;
    private static final double d1 = (1 - Math.exp(-1));
    private static final double f2 = Math.log(Math.sqrt(2*Math.PI));

//    public SpeciesTreeSimplePrior(SpeciesTreeModel sTree, ParametricDistributionModel dist, ParametricDistributionModel tipsPrior) {
//        super(new CompoundModel("STprior"));
//        this.sTree = sTree;
//        this.dist = dist;
//        this.tips = tipsPrior;
//
//        final CompoundModel cm = (CompoundModel)this.getModel();
//        cm.addModel(tipsPrior);
//        cm.addModel(dist);
//        cm.addModel(sTree);
//    }

    public SpeciesTreeSimplePrior(SpeciesTreeModel sTree, Parameter sigma, ParametricDistributionModel tipsPrior) {
           super(new CompoundModel("STprior"));
           this.sTree = sTree;
           this.sigma = sigma;
           this.tips = tipsPrior;

           final CompoundModel cm = (CompoundModel)this.getModel();
           cm.addModel(tipsPrior);
           cm.addModel(sTree);
       }

    protected double calculateLogLikelihood() {
        double ll = 0;
        final NodeRef root = sTree.getRoot();
        double l = sTree.getNodeHeight(root) + ((VDdemographicFunction)sTree.getNodeDemographic(root)).naturalLimit() ;

        for(int nn = 0; nn < sTree.getNodeCount(); ++nn) {
            final NodeRef n = sTree.getNode(nn);
            final DemographicFunction demog = sTree.getNodeDemographic(n);

            if( sTree.isExternal(n) ) {
                ll += tips.logPdf(demog.getDemographic(0));
            }

            final double branch = sTree.isRoot(n) ? ((VDdemographicFunction)demog).naturalLimit()
                    :  sTree.getBranchLength(n);
            final double avg = branch/demog.getIntegral(0, branch);

            final double p = demog.getDemographic(branch);
            final double s = sigma.getParameterValue(0) * Math.sqrt((1 - Math.exp(-branch/l)) / d1);
            //final double vx = (new LogNormalDistribution(-0.5 * s * s, s)).logPdf(p / avg);

            double z = p/avg;
            final double v2 = Math.log(z);
            final double v1 = v2 / s + s / 2;
            final double v = -Math.log(s) - f2 - v2 - v1*v1/2;
            ll += v;
            //ll += dist.logPdf(p/avg);
//
//            double[] pa = {demog.getDemographic(0), demog.getDemographic(branch)};
//            for( double p : pa ) {
//               ll += dist.logPdf(p/avg);
//            }
        }

        return ll;
    }

    protected boolean getLikelihoodKnown() {
		return false;
	}

}
