/*
 * SplitSystem.java
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

package dr.evolution.tree;

import dr.evolution.util.TaxonList;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * data structure for a set of splits 
 *
 * @version $Id: SplitSystem.java,v 1.5 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Korbinian Strimmer
 */
public class SplitSystem
{
	//
	// Public stuff
	//

	/**
	 * @param taxonList  the list of taxa
	 * @param size     number of splits
	 */
	public SplitSystem(TaxonList taxonList, int size)
	{
		this.taxonList = taxonList;
		
		labelCount = taxonList.getTaxonCount();
		splitCount = size;
		
		splits = new boolean[splitCount][labelCount];
	}

	/** get number of splits */
	public int getSplitCount()
	{		
		return splitCount;
	}

	/** get number of labels */
	public int getLabelCount()
	{		
		return labelCount;
	}

	/** get split vector */
	public boolean[][] getSplitVector()
	{		
		return splits;
	}

	/** get split */
	public boolean[] getSplit(int i)
	{		
		return splits[i];
	}


	/** get taxon list */
	public TaxonList getTaxonList() { return taxonList; }

	/**
	  + test whether a split is contained in this split system
	  * (assuming the same leaf order)
	  *
	  * @param split split
	  */
	public boolean hasSplit(boolean[] split)
	{
		for (int i = 0; i < splitCount; i++)
		{
			if (SplitUtils.isSame(split, splits[i])) return true;
		}
			
		return false;
	}


	/** print split system */
	public String toString()
	{
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		
		for (int i = 0; i < labelCount; i++)
		{
			pw.println(taxonList.getTaxon(i));
		}
		pw.println();
		
		
		for (int i = 0; i < splitCount; i++)
		{
			for (int j = 0; j < labelCount; j++)
			{
				if (splits[i][j] == true)
					pw.print('*');
				else
					pw.print('.');
			}
			
			pw.println();
		}

		return sw.toString();
	}

	
	//
	// Private stuff
	//
	
	private int labelCount, splitCount;
	private TaxonList taxonList;
	private boolean[][] splits;
}
