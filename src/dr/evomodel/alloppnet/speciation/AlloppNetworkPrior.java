/*
 * AlloppNetworkPrior.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.speciation;

import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;


/**
 * 
 * Calculates prior likelihood for an allopolyploid network.
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */





public class AlloppNetworkPrior extends Likelihood.Abstract {
	AlloppSpeciesNetworkModel asnm;
	AlloppNetworkPriorModel prior;
    double numHybsLogL[];



    public AlloppNetworkPrior(AlloppNetworkPriorModel prior, AlloppSpeciesNetworkModel asnm) {
		super(prior);
		this.asnm = asnm;
		this.prior = prior;
        asnm.setHybPopModel(prior.getHybridPopModel());

        int ndips = asnm.getNofDiploids();
        numHybsLogL = new double[asnm.maxNumberOfHybPopParameters()+1];
        for (int h = 0; h < numHybsLogL.length; h++) {
            numHybsLogL[h] = -h * Math.log(7.0*Math.sqrt(ndips));
            // grjtodo-soon the form of the function is experimental
        }

		
	}
	
	

	@Override
	protected boolean getLikelihoodKnown() {
		return false;
	}

	@Override
	protected double calculateLogLikelihood() {
        double llhood = 0.0;
        //network topology and times prior
        llhood += loglikelihoodEvents();
        //System.out.print(llhood); System.out.print(" ");
        llhood += loglikeNumHybridizations();
        //System.out.print(llhood); System.out.print(" ");
        // population prior for tips
        Parameter tippvals = asnm.getTipPopValues();
        ParametricDistributionModel tipmodel = prior.getTipPopModel();
        for (int i = 0; i < tippvals.getDimension(); i++) {
            llhood += tipmodel.logPdf(tippvals.getParameterValue(i));
        }
        // population prior for root ends
        Parameter rootpvals = asnm.getRootPopValues();
        ParametricDistributionModel rootmodel = prior.getRootPopModel();
        for (int i = 0; i < rootpvals.getDimension(); i++) {
            llhood += rootmodel.logPdf(rootpvals.getParameterValue(i));
        }
        // population prior for new hybrids
        ParametricDistributionModel hybmodel = prior.getHybridPopModel();
        for (int i = 0; i < asnm.getNumberOfTetraTrees(); i++) {
            llhood += hybmodel.logPdf(asnm.getOneHybPopValue(i));
        }
        //System.out.println(llhood);
        return llhood;
	}





    private double loglikeNumHybridizations() {
        int nhybs = asnm.getNumberOfTetraTrees();
        //System.out.print(nhybs); System.out.print(" ");
        return numHybsLogL[nhybs];
    }


    /*
       * Going backwards in time this gives probabilities to three types
       * of events: diploid-diploid joins, tet-tet joins, and hybridization events.
       */
    private double loglikelihoodEvents() {

        double lambda = prior.getRate().getParameterValue(0);

        ArrayList<Double> heights = new ArrayList<Double>();

        AlloppDiploidHistory adhist = asnm.getDiploidHistory();
        adhist.collectInternalAndHybHeights(heights);
        int nttrees = asnm.getNumberOfTetraTrees();
        for (int tt = 0; tt < nttrees; tt++) {
            AlloppLeggedTree ttree = asnm.getTetraploidTree(tt);
            ttree.collectInternalHeights(heights);
        }
        double loglhood = 0.0;
        for (double t : heights) {
            loglhood += logexpPDF(t, lambda);
        }
        return loglhood;
    }


    private double logexpPDF(double x, double rate) {
        return Math.log(rate) - rate*x;
    }

}
