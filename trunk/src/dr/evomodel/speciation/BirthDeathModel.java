/*
* BirthDeathModel.java
*
* Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;


/**
 * This class implements the Yang and Rannala (1997) birth-death speciation model.
 * @author Andrew Rambaut
 */
public class BirthDeathModel extends SpeciationModel {

	public static final String BIRTH_DEATH_MODEL = "birthDeathModel";
	public static String BIRTH_RATE = "birthRate";
	public static String DEATH_RATE = "deathRate";
	public static String SAMPLING_PROPORTION = "samplingProportion";


	public BirthDeathModel(Parameter birthRateParameter, Parameter deathRateParameter, Parameter samplingProportionParameter, Type units) {

		super(BIRTH_DEATH_MODEL, units);

		this.birthRateParameter = birthRateParameter;
		addParameter(birthRateParameter);
		birthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.deathRateParameter = deathRateParameter;
		addParameter(deathRateParameter);
		deathRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

		this.samplingProportionParameter = samplingProportionParameter;
		addParameter(samplingProportionParameter);
		samplingProportionParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
	}

	public double getBirthRate() { return birthRateParameter.getParameterValue(0); }
	public void setBirthRate(double birthRate) { birthRateParameter.setParameterValue(0, birthRate); }

	public double getDeathRate() { return deathRateParameter.getParameterValue(0); }
	public void setDeathRate(double deathRate) { deathRateParameter.setParameterValue(0, deathRate); }

	public double getSamplingProportion() { return samplingProportionParameter.getParameterValue(0); }
	public void setSamplingProportion(double samplingProportion) { samplingProportionParameter.setParameterValue(0, samplingProportion); }

	//
	// functions that define a speciation model
	//
	public double logTreeProbability(int taxonCount) {
		double lambda = getBirthRate();
		double mu = getDeathRate();
		if (Math.abs(lambda - mu) < 1e-20) {
			return Double.NEGATIVE_INFINITY;
		} else {
			return (taxonCount - 2.0) * Math.log(lambda);
		}
	}

	//
	// functions that define a speciation model
	//
    public double logNodeProbability(Tree tree, NodeRef node) {
        if (tree.getRoot() == node) return 0.0;

        double nodeHeight = tree.getNodeHeight(node);
        double rootHeight = tree.getNodeHeight(tree.getRoot());

		double lambda = getBirthRate();
		double mu = getDeathRate();
		double rho = getSamplingProportion();

        // This is probably better done by including a test prior, if required
        if (lambda < mu) return Double.NEGATIVE_INFINITY;

		return logPDF(nodeHeight, rootHeight, lambda, mu, rho);
	}

	private double logPDF(double nodeTime, double rootTime, double lambda, double mu, double rho)
	{
//   this calculates the kernel density from the birth-death process with
//   species sampling.

		final double small=1e-20;

		if (Math.abs(lambda - mu) < small) {
			double term = 1 + rho * lambda * nodeTime;
			return Math.log((1+rho*lambda * rootTime)/(rootTime * term * term));
		}

		double expmlt1 = Math.exp((mu - lambda) * nodeTime);
		double term1 = p0tBD(lambda, mu, rho, expmlt1);

		double expmlt0 = Math.exp((mu - lambda) * rootTime);
		double term2 = 1.0 - p0tBD(lambda, mu, rho, expmlt0) / rho * expmlt0;

        return Math.log(term1 * term1 * lambda / term2 * rho * expmlt1);

//		return 2.0 * Math.log(term1) + Math.log(lambda) - Math.log(term2 * rho * expmlt1);
	}

	private double p0tBD(double lambda, double mu, double rho, double expmlt) {
		return rho * (lambda - mu) / (rho * lambda + (lambda * (1.0 - rho) - mu) * expmlt );
	}

	/*
		 This is the original code from mcmctree.c, part of the PAML package written
		 by Ziheng Yang.

    #define P0t_BD(expmlt) (rho*(lambda-mu)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt))

	double PDFkernelBD(double t, double t1, double vt1, double lambda, double mu, double rho)
	{
//   this calculates the kernel density from the birth-death process with
//   species sampling.

	   double pdf, expmlt, small=1e-20;

	   if(fabs(lambda-mu)<small)
	      pdf = (1+rho*lambda*t1)/(t1*square(1+rho*lambda*t));
	   else {
	      expmlt=exp((mu-lambda)*t);
	      pdf = P0t_BD(expmlt);
	      pdf = pdf*pdf * lambda/(vt1*rho) * expmlt;
	   }
	   return(pdf);
	}

	double CDFkernelBD(double t, double t1, double vt1, double lambda, double mu, double rho)
	{
//   this calculates the CDF for the kernel density from the birth-death process with
//   species sampling.

	   double cdf, expmlt, small=1e-20;

	   if(fabs(lambda-mu)<small)
	      cdf = (1+rho*lambda*t1)*t/(t1*(1+rho*lambda*t));
	   else {
	      expmlt=exp((mu-lambda)*t);
	      cdf = rho*lambda/vt1 * (1-expmlt)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt);
	   }
	   return(cdf);
	}

	double lnpriorTimes (void)
	{
//   This calculates the prior density of node times in the master species tree:
//   sptree.nodes[].age.  It uses sptree.nodes[].tfossil, sptree.nodes[].fossil[],
//   and lamdba, mu, rho from the birth-death process with species sampling
//   (data.birth, data.death, data.sampling).
//   The routine sorts the node ages in the species tree and then uses the
//   birth-death prior conditional on the calibration points.  t[0] is t1 in the
//   paper.
//
//   rank[1]=3: means that age of node [ns+1] is the 3rd youngest.
//   nodesort[3]=1 means that the 3rd yougest age is node number ns+1.
//   Root (node ns) is excluded in the ranking.

	   int  i,j,k, n1=sptree.nspecies-1, rank[NS-1], rankprev, nfossil;
	   int  nodesort[NS-1]; // nodes with times sorted
	   double t[NS-1], lnpt, Scale, expmlt=1, vt1, P0t1, cdf, cdfprev, small=1e-20;

	   double lambda=data.birth, mu=data.death, rho=data.sampling;

	   if(sptree.root!=sptree.nspecies) error2("node number for root fixed.");
	   for(j=sptree.nspecies; j<sptree.nnode; j++)
	      { t[j-sptree.nspecies]=sptree.nodes[j].age; }

	   // ranking the (n-2) node ages
	   for(i=1,rank[0]=-1,nodesort[0]=-1;i<n1;i++) {
	      for(j=1,k=1;j<n1;j++)  if(j!=i && t[j]<=t[i]) k++;
	      rank[i]=k;
	      nodesort[k]=i;
	   }

	   if(debug==1) {
	      matout2(F0, t, 1, n1, 9, 5);
	      FOR(j,n1) printf("%9d", rank[j]);  FPN(F0);
	      FOR(j,n1) printf("%9d", nodesort[j]);  FPN(F0);
	   }

	   // calculate vt1, needed only if (lambda!=mu)
	   if(fabs(lambda-mu)>small) {
	      expmlt= exp((mu-lambda)*t[0]);
	      P0t1  = P0t_BD(expmlt);
	      vt1   = 1-P0t1/rho*expmlt;
	   }
	   else {
	      P0t1 = rho/(1+rho*mu*t[0]);
	      vt1  = mu*t[0]*P0t1;
	   }
	   // calculate f_BD(t_{_C)}, joint of the remaining times
	   for(j=1,lnpt=1,Scale=0; j<n1; j++) {
	      if(!sptree.nodes[sptree.nspecies+j].fossil)
	         lnpt *= PDFkernelBD(t[j], t[0], vt1, lambda, mu, rho);
	      if(j%50==0) { Scale+=log(lnpt); lnpt=1; }
	   }
	   lnpt=Scale+log(lnpt);

//	      Now calculate f_BD(t_C), marginal for calibration nodes.
//	      This goes through the nodes in the order of their ages, so that node j
//	      is the k-th youngest.

	   for(k=1,nfossil=0,rankprev=0,cdfprev=0; k<n1; k++) {
	      if(!sptree.nodes[j=sptree.nspecies+nodesort[k]].fossil)
	         continue;
	      cdf = CDFkernelBD(t[nodesort[k]], t[0], vt1, lambda, mu, rho);
	      if(k-rankprev-1>0) {
	         lnpt -= (k-rankprev-1)*log(cdf-cdfprev);
	         lnpt += LnGamma((double)k-rankprev-1+1);
	      }
	      if(debug==1)
	         printf("Fossil at node %d age %9.4f  rank diff %d - %d cdf %9.5f\n", j,t[nodesort[k]], k, rankprev, cdf);
	      rankprev=k;  cdfprev=cdf;
	      nfossil++;
	   }
	   if(nfossil && n1-1-rankprev>0) {
	      lnpt -= (n1-1.-rankprev)*log(1-cdfprev);
	      lnpt += LnGamma(n1-1.-rankprev+1);
	   }

	   // Adhockery, added 3 May 2004, in Edmonton
	   if(!sptree.nodes[sptree.root].fossil)
	      lnpt += 2*log(P0t1*(1-vt1))+(n1-1)*log(vt1);

	   lnpt += lnptC_Fossil();
	   if(debug==1) printf("\npdf = %.12f\n", exp(lnpt));

	   return (lnpt);
	}
*/



	/**
	 * Parses an element from an DOM document into a SpeciationModel. Recognises
	 * birthDeathModel.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return BIRTH_DEATH_MODEL; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Type units = XMLParser.Utils.getUnitsAttr(xo);

			Parameter birthParameter = (Parameter)xo.getSocketChild(BIRTH_RATE);
			Parameter deathParameter = (Parameter)xo.getSocketChild(DEATH_RATE);
			Parameter samplingParameter = (Parameter)xo.getSocketChild(SAMPLING_PROPORTION);

			return new BirthDeathModel(birthParameter, deathParameter, samplingParameter, units);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "The Yang & Rannala (1997) model of speciation.";
		}

		public Class getReturnType() { return BirthDeathModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new ElementRule(BIRTH_RATE, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
				new ElementRule(DEATH_RATE, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
				new ElementRule(SAMPLING_PROPORTION, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
				XMLUnits.SYNTAX_RULES[0]
		};
	};


	private final Parameter birthRateParameter;
	private final Parameter deathRateParameter;
	private final Parameter samplingProportionParameter;
}


