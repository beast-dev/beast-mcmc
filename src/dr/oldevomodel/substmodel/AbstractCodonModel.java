/*
 * AbstractCodonModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;

/**
 * base class for codon rate matrices
 *
 * @version $Id: AbstractCodonModel.java,v 1.6 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Korbinian Strimmer
 */
abstract public class AbstractCodonModel extends AbstractSubstitutionModel  {

	// Constructor
	public AbstractCodonModel(String name, Codons codonDataType, FrequencyModel freqModel)
	{
		super(name, codonDataType, freqModel);
		
		this.codonDataType = codonDataType;
		this.geneticCode = codonDataType.getGeneticCode();
	}
	
	protected void frequenciesChanged() {
		// Nothing to precalculate
	}
	
	protected void ratesChanged() {
		// Nothing to precalculate
	}
	
	protected Codons codonDataType;
	protected GeneticCode geneticCode;
}