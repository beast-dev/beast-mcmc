/*
 * ARGPartitionLikelihood.java
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

package dr.evomodel.arg;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;

public abstract class ARGPartitionLikelihood extends AbstractModelLikelihood {

	
	private final ARGModel arg;
	
	public abstract double[] generatePartition();
	public abstract double getLogLikelihood(double[] partition);
	
	public double getLogLikelihood() {
    			
    	double logPrior = 0;
    	
    	for(int i = 0, n = getReassortmentNodeCount(); i < n; i++){
    		logPrior += getLogLikelihood(getPartition(i));
    	}
    	
    	return logPrior;
	}
	
	public double[] getPartition(int i){
		if(arg.getReassortmentNodeCount() == 0){
			return null;
		}
		
		return arg.getPartitioningParameters().getParameter(i).getParameterValues();
	}
	
	public ARGPartitionLikelihood(String id, ARGModel arg){
		super(id);
		
		
		this.arg = arg;
	}
	
	public int getNumberOfPartitionsMinusOne(){
		return arg.getNumberOfPartitions() - 1;
	}
	
	public int getReassortmentNodeCount(){
		return arg.getReassortmentNodeCount();
	}
	
	public Model getModel() {
		return this;
	}
	
	public void makeDirty() {
		
	}
}
