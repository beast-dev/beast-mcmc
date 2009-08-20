/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import java.util.List;

import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.RelativeRatesType;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class SubstitutionModelOptions extends ModelOptions {

	// Instance variables
    private final BeautiOptions options;
   
    

    public SubstitutionModelOptions(BeautiOptions options) {    	
    	this.options = options;
               
        initGlobalSubstitutionModelParaAndOpers();
    }
    
    private void initGlobalSubstitutionModelParaAndOpers() {    	
    	
    	// A vector of relative rates across all partitions...
        createParameter("allMus", "All the relative rates regarding codon positions");

        // This only works if the partitions are of the same size...
//      createOperator("centeredMu", "Relative rates",
//              "Scales codon position rates relative to each other maintaining mean", "allMus",
//              OperatorType.CENTERED_SCALE, 0.75, 3.0);
        createOperator("deltaMu", RelativeRatesType.MU_RELATIVE_RATES.toString(),
        		 "Currently use to scale codon position rates relative to each other maintaining mean", "allMus",
                OperatorType.DELTA_EXCHANGE, 0.75, 3.0);
    	
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {    	    	
	
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
    	if (hasCodon()) {
            Operator deltaMuOperator = getOperator("deltaMu");

            // update delta mu operator weight
            deltaMuOperator.weight = 0.0;
            for (PartitionSubstitutionModel pm : options.getPartitionSubstitutionModels()) {
                deltaMuOperator.weight += pm.getCodonPartitionCount();
            }

            ops.add(deltaMuOperator);
        }
    }
    
    /////////////////////////////////////////////////////////////
    /**
     * @return true either if the options have more than one partition or any partition is
     *         broken into codon positions.
     */
    public boolean hasCodon() {
//        final List<PartitionSubstitutionModel> models = options.getPartitionSubstitutionModels();
//        return (models.size() > 1 || models.get(0).getCodonPartitionCount() > 1);
    	for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            if (model.getCodonPartitionCount() > 1) {
            	return true;
            }
        }
        return false;
    }
    
	@Override
	public String getPrefix() {
		// TODO Auto-generated method stub
		return null;
	}


}
