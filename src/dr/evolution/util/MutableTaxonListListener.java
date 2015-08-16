/*
 * MutableTaxonListListener.java
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

package dr.evolution.util;

/**
 * Interface for a listener to a mutable list of taxa.
 *
 * @version $Id: MutableTaxonListListener.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public interface MutableTaxonListListener {

	/** Called when a taxon has been added to the taxonList */
	void taxonAdded(TaxonList taxonList, Taxon taxon);
	
	/** Called when a taxon has been successfully removed from the taxonList */
	void taxonRemoved(TaxonList taxonList, Taxon taxon);
	
	/** Called when one or more taxon has been edited */
	void taxaChanged(TaxonList taxonList);
}
