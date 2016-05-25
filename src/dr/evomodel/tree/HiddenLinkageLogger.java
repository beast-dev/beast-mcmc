/*
 * HiddenLinkageLogger.java
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

package dr.evomodel.tree;

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

/**
 * @author Aaron Darling
 */
public class HiddenLinkageLogger extends MCLogger {

	HiddenLinkageModel hlm;
	public HiddenLinkageLogger(HiddenLinkageModel hlm, LogFormatter formatter, int logEvery) 
	{
		super(formatter, logEvery, false);
		this.hlm = hlm;
	}

    public void startLogging() {
    	String[] labels = new String[1 + hlm.getData().getReadsTaxa().getTaxonCount()];
    	labels[0] = "iter";
    	for(int i=1; i<labels.length; i++){
    		labels[i] = hlm.getData().getReadsTaxa().getTaxonId(i-1);
    	}
    	this.logLabels(labels);
    }
    public void log(long state) {
    	if(state % logEvery != 0)
    		return;
    	String[] values = new String[1 + hlm.getData().getReadsTaxa().getTaxonCount()];
    	values[0] = new Long(state).toString();
    	for(int i=1; i<values.length; i++){
    		values[i] = new Integer(hlm.getLinkageGroupId(hlm.getData().getReadsTaxa().getTaxon(i-1))).toString();
    	}
    	this.logValues(values);
    }
    
}
