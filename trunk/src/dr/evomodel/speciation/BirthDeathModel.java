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

import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;


/**
 * This class implements the Yang and Rannala (1997) birth-death speciation model.
 *@author Andrew Rambaut
 */
public class BirthDeathModel extends SpeciationModel {

    public static final String BIRTH_DEATH_MODEL = "birthDeathModel";
    public static String BIRTH_RATE = "birthRate";
    public static String DEATH_RATE = "deathRate";
    public static String SAMPLING_PROPORTION = "samplingProportion";


    public BirthDeathModel(Parameter birthRateParameter, Parameter deathRateParameter, Parameter samplingProportionParameter, int units) {

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
    public double logNodeHeightProbability(double nodeHeight, double rootHeight) {
        double lambda = getBirthRate();
        double mu = getDeathRate();
        double rho = getSamplingProportion();

        if (lambda < mu) {
            return Double.NEGATIVE_INFINITY;
        }

        // It looks as if Yang and Rannala (1997) have the root time constrained to be 1.0. 
        double t1 = nodeHeight / rootHeight;
        double t0 = 1.0;
//        double t1 = nodeHeight;
//        double t0 = rootHeight;

        double lnP = 0.0;

        if (Math.abs(lambda - mu) < 1e-20) {
            lnP = Math.log (1.0 + rho * mu) - (2.0 * Math.log(1.0 + rho * mu * t1));
        } else {
            double p0t1 = rho * (lambda - mu) / (rho * lambda + (lambda * (1.0 - rho) - mu) * Math.exp((mu - lambda) * t1) );
            lnP = (Math.log(1.0 / rho) + 2.0 * Math.log(p0t1) + (mu - lambda) * t1);

            double p0t0 =  rho *(lambda - mu) / (rho * lambda + (lambda * (1.0 - rho) - mu) * Math.exp((mu - lambda) * t0) );
            lnP = lnP - (Math.log(1.0 - (1.0 / rho) * p0t0 * Math.exp((mu - lambda) * t0)));
        }

        return lnP;
    }

    /*
     This is the original code from mcmc.c, part of the MrBayes program written
     by John P. Huelsenbeck & Fredrik Ronquist. MrBayes is distributed under a
     GNU GPL license - see original software for details.
     
    int LnBirthDeathPriorPr (Tree *t, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF)
    {
        int				i, j, nNodes;
        MrBFlt			rootTime=0.0, *nt;
        TreeNode		*p;

        nt = (MrBFlt *)SafeMalloc((size_t) (t->nIntNodes) * sizeof(MrBFlt));
        if (!nt)
            {
            printf ("\n   ERROR: Problem allocating nt\n");
            return (ERROR);
            }

        // get the node times and put them into a vector
        for (i=j=0; i<t->nIntNodes; i++)
            {
            p = t->intDownPass[i];
            if (p->anc->anc != NULL)
                nt[j++] = p->nodeDepth;
            else
                rootTime = p->nodeDepth;
            }
        nNodes = j;

        // rescale all of the node times on the tree
        for (i=0; i<nNodes; i++)
            nt[i] /= rootTime;

        // I think this is correct. It looks as if Yang and Rannala (1997)
           have the root time constrained to be 1.0.
        rootTime = 1.0;

        / calculate probabilities of tree
        if (AreDoublesEqual(sR,eR,ETA)==NO)
            {
            (*prob) = (numLocalTaxa - 2.0) * log(sR);
            for (i=0; i<nNodes; i++)
                (*prob) += LnP1 (nt[i], sR, eR, sF) - LnVt (rootTime, sR, eR, sF);
            }
        else
            {
            (*prob) = 0.0;
            for (i=0; i<nNodes; i++)
                (*prob) += log (1.0 + sF * eR) - (2.0 * log(1.0 + sF * eR * nt[i]));
            }

        free (nt);

        return (NO_ERROR);
    }

    MrBFlt LnP1 (MrBFlt t, MrBFlt l, MrBFlt m, MrBFlt r)
    {

        MrBFlt		p0t;
        p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
        return (log(1.0/r) + 2.0*log(p0t) + (m-l)*t);

    }

    MrBFlt LnVt (MrBFlt t, MrBFlt l, MrBFlt m, MrBFlt r)
    {
        MrBFlt		p0t;
        p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
        return (log(1.0 - (1.0/r) * p0t * exp((m-l)*t)));

    }
    */

    /*
     This is the original code from mcmctree.c, part of the PAML package written
     by Ziheng Yang.

     // calculate vt1, needed only if (lambda!=mu)
     if(fabs(lambda-mu)>small) {
        expmlt=exp((mu-lambda)*t[0]);
        vt1 = 1-P0t(expmlt)/rho*expmlt;
     }
     // calculate f_BD(t_{_C)}, joint of the remaining times
     for(j=1,lnpt=1,Scale=0; j<n1; j++) {
        if(!sptree.nodes[sptree.nspecies+j].fossil)
           lnpt *= PDFkernelBD(t[j], t[0], vt1, lambda, mu, rho);
        if(j%50==0) { Scale+=log(lnpt); lnpt=1; }
     }

    #define P0t(expmlt) (rho*(lambda-mu)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt))

    double PDFkernelBD(double t, double t1, double vt1, double lambda, double mu, double rho)
    {
        // this calculates the kernel density from the birth-death process with
        //species sampling.
       double pdf, expmlt, small=1e-20;

       if(fabs(lambda-mu)<small)
          pdf = (1+rho*lambda*t1)/(t1*square(1+rho*lambda*t));
       else {
          expmlt=exp((mu-lambda)*t);
          pdf = P0t(expmlt);
          pdf = pdf*pdf * lambda/(vt1*rho) * expmlt;
       }
       return(pdf);
    }

    double CDFkernelBD(double t, double t1, double vt1, double lambda, double mu, double rho)
    {
        // this calculates the CDF for the kernel density from the birth-death process with
        // species sampling.
       double cdf, expmlt, small=1e-20;

       if(fabs(lambda-mu)<small)
          cdf = (1+rho*lambda*t1)*t/(t1*(1+rho*lambda*t));
       else {
          expmlt=exp((mu-lambda)*t);
          cdf = rho*lambda/vt1 * (1-expmlt)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt);
       }
       return(cdf);
    }
    */

    /**
     * Parses an element from an DOM document into a SpeciationModel. Recognises
     * birthDeathModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return BIRTH_DEATH_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int units = XMLParser.Utils.getUnitsAttr(xo);

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
