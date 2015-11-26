/*
 * ARGDistinctTreeCountStatistic.java
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

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class ARGDistinctTreeCountStatistic extends Statistic.Abstract{

	public static final String ARG_TREE_COUNT = "argTreeCount";
	public static final String PRINT_TREES = "printTrees";
	
	private ARGModel arg;
	private int numberOfPartitions;
	
	public ARGDistinctTreeCountStatistic(ARGModel arg){
		this.arg = arg;
		numberOfPartitions = arg.getNumberOfPartitions();
	}
	
	public int getDimension() {
		
		return 1;
	}

	
	public double getStatisticValue(int dim){
		ArrayList<String> listOfTrees = new ArrayList<String>(numberOfPartitions/10);
		
		for(int i = 0; i < numberOfPartitions; i++){
			ARGTree tree = new ARGTree(arg,i);
			
			String newick = tree.getUniqueNewick();
			
			
			
			if(!listOfTrees.contains(newick)){
				listOfTrees.add(newick);
				
			}
		}
		
		
		
			
		return listOfTrees.size();
	}
	
	private class DistinctTreeCount implements Comparable<DistinctTreeCount>{
		private final String tree;
		private int count;
		
		public DistinctTreeCount(String tree, int count){
			this.tree = tree;
			this.count = count;
		}
		
		public void increaseCount(){
			count++;
		}

		public int compareTo(DistinctTreeCount o) {
			if(this.tree.equals(o.tree)){
				return 0;
			}else if(this.count < o.count){
				return 1;
			}
			return -1;
		}
		
		public boolean equals(Object o){
			try{
				DistinctTreeCount dtc = (DistinctTreeCount)o;
				
				if(dtc.compareTo(this) == 0){
					return true;
				}
				return false;
			}catch(ClassCastException e){
				return false;
			}
			
		}
		
		public String toString(){
			return tree + " " + count;
		}
		
		
		
		
	};
	
	
	
	public String getFullOutput(){
		ArrayList<String> listOfTrees = new ArrayList<String>(numberOfPartitions);
		ArrayList<Integer> numbers = new ArrayList<Integer>(numberOfPartitions);
		
		for(int i = 0; i < numberOfPartitions; i++){
			ARGTree tree = new ARGTree(arg,i);
			
			String newick = tree.getUniqueNewick();
			if(!listOfTrees.contains(newick)){
				listOfTrees.add(newick);
				numbers.add(1);
			}else{
				int index = listOfTrees.indexOf(newick);
				numbers.set(index, numbers.get(index) + 1);
			}
		}
		
		ArrayList<DistinctTreeCount> trees = new ArrayList<DistinctTreeCount>(numberOfPartitions);
		
		for(int i = 0; i < listOfTrees.size(); i++){
			trees.add(new DistinctTreeCount(listOfTrees.get(i),numbers.get(i)));
		}
		
		Collections.sort(trees);
		
		return trees.toString();
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
		
			return "Provides number of distinct loci trees";
		}

		public Class getReturnType() {
			return ARGDistinctTreeCountStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
					new ElementRule(ARGModel.class),
			};
		}
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);
			
			return new ARGDistinctTreeCountStatistic(arg);
		}

		public String getParserName() {
			return ARG_TREE_COUNT;
		}
		
	};

	
	
}
