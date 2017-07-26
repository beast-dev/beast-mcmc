/*
 * DriftedTreeClusterLocationsStatistic.java
 *
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

package dr.evomodel.antigenic.phyloclustering.statistics;

import dr.inference.model.*;
import dr.xml.*;

/**
 * @author Charles Cheung
 * 
 * adapted from the code:
 * @author Trevor Bedford
 * @author Marc A. Suchard
 */
public class DriftedTreeClusterLocationsStatistic extends Statistic.Abstract implements VariableListener {

    public static final String DRIFTED_TREE_CLUSTER_LOCATIONS_STATISTIC = "driftedTreeClusterLocationsStatistic";

    public DriftedTreeClusterLocationsStatistic(MatrixParameter locationsParameter, Parameter mu1Scale, Parameter mu2Scale) {
    //public DriftedTreeClusterLocationsStatistic(MatrixParameter locationsParameter, Parameter locationDriftParameter, Parameter mu1Scale, Parameter mu2Scale) {
        this.locationsParameter = locationsParameter;
        locationsParameter.addParameterListener(this);
    //    this.locationDriftParameter = locationDriftParameter;
     //   locationDriftParameter.addParameterListener(this);
        
        this.mu1ScaleParameter = mu1Scale;
        mu1ScaleParameter.addParameterListener(this);
        this.mu2ScaleParameter = mu2Scale;
        mu2ScaleParameter.addParameterListener(this);
    }

    public int getDimension() {
        return locationsParameter.getDimension();
    }

    // strain index
    public int getColumnIndex(int dim) {
        return dim / locationsParameter.getRowDimension();
    }

    // dimension index
    public int getRowIndex(int dim) {
        int x = getColumnIndex(dim);
        return dim - x * locationsParameter.getRowDimension();
    }

    public double getStatisticValue(int dim) {

        double val;

        // x is location count, y is location dimension
        int x = getColumnIndex(dim);
        int y = getRowIndex(dim);
        Parameter loc = locationsParameter.getParameter(x);

       	if (y == 0) {
            //val = loc.getParameterValue(y)* locationDriftParameter.getParameterValue(0);
        		val = loc.getParameterValue(y)* mu1ScaleParameter.getParameterValue(0);
        	}
        else {
            //val = loc.getParameterValue(y);
        	val = loc.getParameterValue(y)* mu2ScaleParameter.getParameterValue(0);
        }
        

        return val;

    }

    public String getDimensionName(int dim) {
        return locationsParameter.getDimensionName(dim);
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String LOCATIONS = "locations";
        //public final static String LOCATION_DRIFT = "locationDrift";

        public final static String MU1_SCALE_PARAMETER = "mu1Scale";
        public final static String MU2_SCALE_PARAMETER = "mu2Scale";
        
        public String getParserName() {
            return DRIFTED_TREE_CLUSTER_LOCATIONS_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter locations = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);
            //Parameter locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
            Parameter mu1Scale = null;
            if (xo.hasChildNamed(MU1_SCALE_PARAMETER)) {
            	mu1Scale = (Parameter) xo.getElementFirstChild(MU1_SCALE_PARAMETER);
            }
            
            Parameter mu2Scale = null;
            if (xo.hasChildNamed(MU2_SCALE_PARAMETER)) {
            	mu2Scale = (Parameter) xo.getElementFirstChild(MU2_SCALE_PARAMETER);
            }  
            
  //          return new DriftedTreeClusterLocationsStatistic(locations, locationDrift, mu1Scale, mu2Scale);
            return new DriftedTreeClusterLocationsStatistic(locations, mu1Scale, mu2Scale);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return DriftedTreeClusterLocationsStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LOCATIONS, MatrixParameter.class),
          //  new ElementRule(LOCATION_DRIFT, Parameter.class, "TO BECOME OBSOLETE", true),
            new ElementRule(MU1_SCALE_PARAMETER, Parameter.class, "Optional parameter for scaling the first dimension of mu"),
            new ElementRule(MU2_SCALE_PARAMETER, Parameter.class, "Optional parameter for scaling the second dimension of mu"),        };
    };

    private MatrixParameter locationsParameter;
//    private Parameter locationDriftParameter;
    private Parameter mu1ScaleParameter;
    private Parameter mu2ScaleParameter;
}
