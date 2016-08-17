/*
 * PathLikelihood.java
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

package dr.inference.model;

import dr.evomodel.continuous.SoftThresholdLikelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * A likelihood function which is simply the product of a set of likelihood functions.
 *
 * @author Andrew Rambaut
 * @author Alex Alekseyenko
 * @author Marc Suchard
 * @version $Id: CompoundLikelihood.java,v 1.19 2005/05/25 09:14:36 rambaut Exp $
 */
public class PathLikelihood implements Likelihood {

    public static final String PATH_LIKELIHOOD = "pathLikelihood";
    public static final String PATH_PARAMETER = "theta";
    public static final String DIFFERENCE = "delta";
    public static final String SOURCE = "source";
    public static final String DESTINATION = "destination";
    public static final String PSUEDO_SOURCE = "sourcePseudoPrior";
    public static final String PSUEDO_DESTINATION = "destinationPseudoPrior";

    public PathLikelihood(Likelihood source, Likelihood destination) {
        this(source, destination, null, null);
    }

    public PathLikelihood(Likelihood source, Likelihood destination,
                          Likelihood pseudoSource, Likelihood pseudoDestination) {
        this.source = source;
        this.destination = destination;
        this.pseudoSource = pseudoSource;
        this.pseudoDestination = pseudoDestination;
        if(source!=null) {
            for (int j = 0; j < source.getModel().getModelCount(); j++) {
                for (int i = 0; i < source.getModel().getModel(j).getModelCount(); i++) {
                    if (source.getModel().getModel(j).getModel(i) instanceof SoftThresholdLikelihood) {
                        thresholdSofteners.add((SoftThresholdLikelihood) source.getModel().getModel(j).getModel(i));
                    }
                }
            }
        }

        if(pseudoSource!=null){
            for (int i = 0; i < pseudoSource.getModel().getModelCount(); i++) {
                if(pseudoSource.getModel().getModel(i) instanceof SoftThresholdLikelihood){
                    thresholdSofteners.add((SoftThresholdLikelihood) pseudoSource.getModel().getModel(i));
                }
            }
        }

        compoundModel.addModel(source.getModel());
        compoundModel.addModel(destination.getModel());
    }

    public double getPathParameter() {
        return pathParameter;
    }

    public void setPathParameter(double pathParameter) {
        this.pathParameter = pathParameter;
        if(thresholdSofteners!=null){
            for(SoftThresholdLikelihood threshold:thresholdSofteners){
                threshold.setPathParameter(pathParameter);
            }
        }
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public Model getModel() {
        return compoundModel;
    }

    public double getLogLikelihood() {
        // Depends on complete model (include pseudo-densities)
        return (source.getLogLikelihood() * pathParameter) + (destination.getLogLikelihood() * (1.0 - pathParameter));
    }

    public Likelihood getSourceLikelihood() {
        return source;
    }

    public Likelihood getDestinationLikelihood() {
        return destination;
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        Set<Likelihood> set = new HashSet<Likelihood>();
        set.add(source);
        set.add(destination);
        return set;
    }

    public void makeDirty() {
        source.makeDirty();
        destination.makeDirty();
    }

    public boolean evaluateEarly() {
        return false;
    }

    public String toString() {

        return Double.toString(getLogLikelihood());

    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new DeltaLogLikelihoodColumn(getId() + "." + DIFFERENCE),
                new PosteriorColumn(getId() + "." + SOURCE),
                new PriorColumn(getId() + "." + DESTINATION),
                new ThetaColumn(getId() + "." + PATH_PARAMETER),
                new PathLikelihoodColumn(getId() + "." + PATH_LIKELIHOOD)
        };
    }
    
    private class PriorColumn extends dr.inference.loggers.NumberColumn {
    	
    	public PriorColumn(String label) {
    		super(label);
    	}
    	
    	public double getDoubleValue() {
    		//assume that the prior is the destination in the BEAST xml
    		return destination.getLogLikelihood();
    	}
    	
    }
    
    private class PosteriorColumn extends dr.inference.loggers.NumberColumn {
    	
    	public PosteriorColumn(String label) {
    		super(label);
    	}
    	
    	public double getDoubleValue() {
    		//assume that the posterior is the source in the BEAST xml
    		return source.getLogLikelihood();
    	}
    	
    }

    private class DeltaLogLikelihoodColumn extends dr.inference.loggers.NumberColumn {

        public DeltaLogLikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            // Remove pseudo-densities
            double totalLikelihoodCorrection = 0;
            if(thresholdSofteners!=null){
                for(SoftThresholdLikelihood threshold:thresholdSofteners){
                    totalLikelihoodCorrection += threshold.getLikelihoodCorrection();
                }
            }

            double logDensity = source.getLogLikelihood() - destination.getLogLikelihood() + totalLikelihoodCorrection;
            if (pseudoSource != null) {
                logDensity -= pseudoSource.getLogLikelihood();
            }
            if (pseudoDestination != null) {
//                System.err.println("value = "+pseudoDestination.getLogLikelihood());
//                logDensity += pseudoDestination.getLogLikelihood();
            }
            return logDensity;
        }
    }

    private class ThetaColumn extends dr.inference.loggers.NumberColumn {

        public ThetaColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return pathParameter;
        }

    }

    private class PathLikelihoodColumn extends dr.inference.loggers.NumberColumn {

        public PathLikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }

    }

    public boolean isUsed() {
        return isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

    public void setId(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    public String prettyName() {
        return Abstract.getPrettyName(this);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PATH_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        	
            Likelihood source = (Likelihood) xo.getElementFirstChild(SOURCE);
            Likelihood destination = (Likelihood) xo.getElementFirstChild(DESTINATION);
            
            Likelihood pseudoSource = null;
            if (xo.hasChildNamed(PSUEDO_SOURCE)) {
                pseudoSource = (Likelihood) xo.getElementFirstChild(PSUEDO_SOURCE);
            }
            Likelihood pseudoDestination = null;
            if (xo.hasChildNamed(PSUEDO_DESTINATION)) {
                pseudoDestination = (Likelihood) xo.getElementFirstChild(PSUEDO_DESTINATION);
            }

            return new PathLikelihood(source, destination, pseudoSource, pseudoDestination);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A likelihood function used for estimating marginal likelihoods and Bayes factors using path sampling.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(SOURCE,
                        new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}),
                new ElementRule(DESTINATION,
                        new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}),
                new ElementRule(PSUEDO_SOURCE,
                        new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}, true),
                new ElementRule(PSUEDO_DESTINATION,
                        new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}, true),
        };

        public Class getReturnType() {
            return CompoundLikelihood.class;
        }
    };

    private final Likelihood source;
    private final Likelihood destination;
    private final Likelihood pseudoSource;
    private final Likelihood pseudoDestination;

    private double pathParameter;

    private final CompoundModel compoundModel = new CompoundModel("compoundModel");

    private ArrayList<SoftThresholdLikelihood> thresholdSofteners=new ArrayList<SoftThresholdLikelihood>();
}