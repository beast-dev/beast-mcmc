/*
 * ARGReassortmentTimingStatistic.java
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Logger;

import dr.evomodel.arg.ARGModel.Node;
import dr.evomodel.arg.operators.ARGPartitioningOperator.PartitionChangedEvent;

import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ARGReassortmentTimingStatistic extends Statistic.Abstract{

	private int dimension;
	private ARGModel arg;
	
	public static final String ARG_TIMING_STATISTIC = "argTimingStatistic";
	public static final String NUMBER_OF_REASSORTMENTS = "reassortments";  //TODO This is probably somewhere else in BEAST.
	
	public ARGReassortmentTimingStatistic(String name, ARGModel arg){
		super(name);
		
		
		this.dimension = arg.getExternalNodeCount() + 1;
		this.arg = arg;
	}
	
	public int getDimension() {
		return dimension;
	}
	
	public String getDimensionName(int dim){
		if(dim == 0){
			return "Root";
		}else if(dim == 1){
			return "RootChild";
		}else if(dim == 2){
			return "ENParent";
		}else if(dim == 3){
			return "DNParent";
		}else if(dim == 4){
			return "DCParent";
		}else if(dim == 5){
			return "CCParent";
		}else if(dim == 6){
			return "FNParent";
		}else if(dim == 7){
			return "ReassortHeight";
		}else if(dim == 8){
			return "CNParent";
		}else if(dim == 9){
			return "CCParentParent";
		}
		
		return "";
	}
	
	public double getStatisticValue(int dim) {
		String max = "((((((<(FC,FN)>,CN),CC),<(FC,FN)>),DC),((EC,EN),DN)),AN);";
		
		
		if(!arg.toExtendedNewick().equals(max)){
			return Double.NaN;
		}
		
		if(dim == 0){
			return arg.getRootHeightParameter().getParameterValue(0);
		}else if(dim == 1){
			Node a = (Node)arg.getRoot();
			
			Node aLeft = a.leftChild;
			Node aRight = a.rightChild;
			
			return Math.max(aLeft.heightParameter.getParameterValue(0), aRight.heightParameter.getParameterValue(0));
		}else if(dim == 2){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("EN")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.heightParameter.getParameterValue(0);
		}else if(dim == 3){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("DN")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.heightParameter.getParameterValue(0);
		}else if(dim == 4){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("DC")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.heightParameter.getParameterValue(0);
		}else if(dim == 5){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("CC")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.heightParameter.getParameterValue(0);
		}else if(dim == 6){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("FN")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.heightParameter.getParameterValue(0);
		}else if(dim == 7){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("FN")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.leftParent.heightParameter.getParameterValue(0);
		}else if(dim == 8){
			int value = 0;
			Node a = (Node)arg.getExternalNode(0);
			
			while(!a.taxon.toString().equals("CN")){
				value++;
				a = (Node)arg.getExternalNode(value);
			}
			return a.leftParent.heightParameter.getParameterValue(0);
		}
		
		int value = 0;
		Node a = (Node)arg.getExternalNode(0);
		
		while(!a.taxon.toString().equals("CC")){
			value++;
			a = (Node)arg.getExternalNode(value);
		}
		return a.leftParent.leftParent.heightParameter.getParameterValue(0);
		
		
		
		
	}
	
	
    
	
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			
			return "";
		}

		public Class getReturnType() {
			
			return ARGReassortmentTimingStatistic.class;
		}

		
		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				new ElementRule(ARGModel.class,false),
				AttributeRule.newStringRule(NAME,true),
			};
		}

		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String name = xo.getId();
//			int dim = xo.getIntegerAttribute(DIMENSION);
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			Logger.getLogger("dr.evomodel").info("Creating timing statistic");
			
			return new ARGReassortmentTimingStatistic(name,arg);
		}

		public String getParserName() {
			return ARG_TIMING_STATISTIC;
		}
		
	};

	
	
}

