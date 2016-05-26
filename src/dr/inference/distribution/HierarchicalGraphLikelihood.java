/*
 * HierarchicalGraphLikelihood.java
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

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Variable;
import dr.math.Binomial;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 *
 *
 *   @author Gabriela Cybis
 */

public class HierarchicalGraphLikelihood extends AbstractModelLikelihood {

    public static final String HIERARCHICAL_GRAPH_LIKELIHOOD = "hierarchicalGraphLikelihood";

    public HierarchicalGraphLikelihood(Parameter hierarchicalIndicator, MatrixParameter strataIndicatorMatrix, Parameter prob) {

        super(HIERARCHICAL_GRAPH_LIKELIHOOD);

        this.hierarchicalIndicator = hierarchicalIndicator;
        this.strataIndicatorMatrix = strataIndicatorMatrix;
        this.prob = prob;
        addVariable(hierarchicalIndicator);
        addVariable(strataIndicatorMatrix);
        addVariable(prob);
        
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************



    public Parameter getHierarchicalIndicator() {
        return this.hierarchicalIndicator;
    }
    

    public MatrixParameter getStrataMatrix() {
        return this.strataIndicatorMatrix;
    }
    
    public Parameter getProb() {
        return this.prob;
    }
    
    public Model getModel() {
        return this;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {

        double p = prob.getParameterValue(0);
        if (p <= 0 || p >= 1) return Double.NEGATIVE_INFINITY;

        double logP = Math.log(p);
        double log1MinusP = Math.log(1.0 - p);

        if ( hierarchicalIndicator.getDimension()!= strataIndicatorMatrix.getRowDimension()) return Double.NEGATIVE_INFINITY;

        double logL = 0.0;
        
        
        for (int j =0; j < strataIndicatorMatrix.getColumnDimension();j++){
       int diff = 0;
        	for (int i = 0; i < hierarchicalIndicator.getDimension(); i++) {
            diff += (int) Math.abs(Math.round(hierarchicalIndicator.getParameterValue(i)-strataIndicatorMatrix.getParameterValue(i,j))); 
        
        } logL += geometricLogLikelihood( diff, logP, log1MinusP);
        
        
        /**   double logL += binomialLogLikelihood(hierarchicalIndicator.getDimension(), diff, logP, log1MinusP);
         *binomialLogLikelihood(hierarchicalIndicator.getDimension(), diff, logP, log1MinusP);
        */
        }
       
        
       
       
     
      
       
       
       
       
   
        return logL;
    }

    
    
    public void makeDirty() {
    }

    public void acceptState() {
        // DO NOTHING
    }

    public void restoreState() {
        // DO NOTHING
    }

    public void storeState() {
        // DO NOTHING
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // DO NOTHING
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // DO NOTHING
    }

    /**
     * @return the bernoulli loglikelihood
     *         when the log of the probability is logP.
     */
    private double binomialLogLikelihood(int trials, int count, double logP, double log1MinusP) {
        return Math.log(Binomial.choose(trials, count)) + (logP * count) + (log1MinusP * (trials - count));
    }

    
    /**
     * @return the geometric loglikelihood
     *         when the log of the probability is logP.
     */
    private double geometricLogLikelihood( int count, double logP, double log1MinusP) {
        return  (log1MinusP ) + (logP * count);
    }
    
    
    
    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

   // Binomial binom = new Binomial();
    Parameter hierarchicalIndicator;
    MatrixParameter strataIndicatorMatrix;
    Parameter prob;
    }





