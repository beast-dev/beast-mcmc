/*
 * ARGTotalLengthStatistic.java
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
import dr.inference.model.Statistic;
import dr.xml.*;

public class ARGTotalLengthStatistic extends Statistic.Abstract{

	public static final String ARG_TOTAL_LENGTH = "argTotalLengthStatistic";
	private ARGModel arg;
	
	public ARGTotalLengthStatistic(String id, ARGModel arg){
		super(id);
		
		this.arg = arg;
	}
	
	public int getDimension() {
		return 1;
	}

	public double getStatisticValue(int dim) {
		
		double length = 0;
		
		for(int i = 0; i < arg.getNodeCount(); i++){
			Node x = (Node) arg.getNode(i);
			
			if(!x.isRoot()){
				length += x.getParent(ARGModel.LEFT).getHeight() - x.getHeight();
				if(x.isReassortment()){
					length += x.getParent(ARGModel.RIGHT).getHeight() - x.getHeight();
				}
			}
		}
		
		return length;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		@Override
		public String getParserDescription() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public Class getReturnType() {
			
			return ARGTotalLengthStatistic.class;
		}

		@Override
		public XMLSyntaxRule[] getSyntaxRules() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String id = xo.getId();
			
			ARGModel a = (ARGModel) xo.getChild(ARGModel.class);
			
			return new ARGTotalLengthStatistic(id,a);
		}

		public String getParserName() {
			// TODO Auto-generated method stub
			return ARG_TOTAL_LENGTH;
		}
		
	};

}
