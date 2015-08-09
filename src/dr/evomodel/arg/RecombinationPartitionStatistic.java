/*
 * RecombinationPartitionStatistic.java
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

import dr.evomodel.arg.ARGModel.Node;
import dr.evomodel.arg.operators.ARGAddRemoveEventOperator;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class RecombinationPartitionStatistic extends Statistic.Abstract{

	public final static String RECOMBINATION_PARTITION_STATISTIC = "partitionStatistic";
	private int dimension;
	private ARGModel arg;
	private String[] taxaNames;
	
	public RecombinationPartitionStatistic(String id, ARGModel arg){
		
		setId(id);
		
		this.arg = arg;
		
		this.dimension = arg.getExternalNodeCount();
		taxaNames = new String[this.dimension];
		for(int i = 0; i < taxaNames.length; i++){
			taxaNames[i] = "" + ((Node)arg.getExternalNode(i)).taxon;
		}
	}
	
	
	public int getDimension() {
		return dimension;
	}

	public String getDimensionName(int dim){
		return "Taxa" + taxaNames[dim];
	}
	
	public double getStatisticValue(int dim) {
		
		Node x = (Node)arg.getExternalNode(dim);
		
		assert x.taxon.toString().equals(taxaNames[dim]);
		
		boolean c = x.hasReassortmentAncestor();
		
		if(c){
			return 1.0;
		}
		return 0.0;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			return RecombinationPartitionStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String id = xo.getId();
			
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			return new RecombinationPartitionStatistic(id,arg);
		}

		public String getParserName() {
			return RECOMBINATION_PARTITION_STATISTIC;
		}
		
	};

}
