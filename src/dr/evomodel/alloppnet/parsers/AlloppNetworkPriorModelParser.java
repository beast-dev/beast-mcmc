/*
 * AlloppNetworkPriorModelParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.parsers;

import dr.evolution.util.Units;
import dr.evomodel.alloppnet.speciation.AlloppNetworkPriorModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.inference.distribution.ParametricDistributionModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */

/*
<alloppNetworkPriorModel id="network.prior.model" units="substitutions" >
        <eventRate>
        <parameter id="apspnetwork.prior.eventRate" value="100" lower="0" upper="Infinity" />
        </eventRate>
        <populationScalingFactor>
        <parameter id="population.scaling.factor" value="0.2" lower="0" upper="Infinity" />
        </populationScalingFactor>
        <tipPopulationDistribution>
        <gammaDistributionModel>
        <shape>
        4
        </shape>
        <scale>
        <parameter idref="population.scaling.factor" />
        </scale>
        </gammaDistributionModel>
        </tipPopulationDistribution>
        <rootPopulationDistribution>
        <gammaDistributionModel>
        <shape>
        2
        </shape>
        <scale>
        <parameter idref="population.scaling.factor" />
        </scale>
        </gammaDistributionModel>
        </rootPopulationDistribution>
        <hybridPopulationDistribution>
        <gammaDistributionModel>
        <shape>
        1
        </shape>
        <scale>
        <parameter idref="population.scaling.factor" />
        </scale>
        </gammaDistributionModel>
        </hybridPopulationDistribution>
        </alloppNetworkPriorModel>
*/
public class AlloppNetworkPriorModelParser extends AbstractXMLObjectParser {

	public static final String ALLOPPNETWORKPRIORMODEL = "alloppNetworkPriorModel";
	public static final String EVENTRATE = "eventRate";
    public static final String POPULATION_SCALING_FACTOR = "populationScalingFactor";
    public static final String TIP_POPULATION_DISTRIBUTION = "tipPopulationDistribution";
    public static final String ROOT_POPULATION_DISTRIBUTION = "rootPopulationDistribution";
    public static final String HYBRID_POPULATION_DISTRIBUTION = "hybridPopulationDistribution";


	public String getParserName() {
		return ALLOPPNETWORKPRIORMODEL;
	}



    @Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
		final XMLObject erXo = xo.getChild(EVENTRATE);
		final Parameter eventrate = (Parameter) erXo.getChild(Parameter.class);
        final XMLObject psfXo = xo.getChild(POPULATION_SCALING_FACTOR);
        final Parameter popscalingfactor = (Parameter) psfXo.getChild(Parameter.class);
        final XMLObject tpdXo = xo.getChild(TIP_POPULATION_DISTRIBUTION);
        ParametricDistributionModel tippopmodel = (ParametricDistributionModel) tpdXo.getChild(ParametricDistributionModel.class);
        final XMLObject rpdXo = xo.getChild(ROOT_POPULATION_DISTRIBUTION);
        ParametricDistributionModel rootpopmodel = (ParametricDistributionModel) rpdXo.getChild(ParametricDistributionModel.class);
        final XMLObject hpdXo = xo.getChild(HYBRID_POPULATION_DISTRIBUTION);
        ParametricDistributionModel hybpopmodel = (ParametricDistributionModel) hpdXo.getChild(ParametricDistributionModel.class);
        return new AlloppNetworkPriorModel(eventrate, popscalingfactor, tippopmodel, rootpopmodel, hybpopmodel, units);
	}



	private XMLSyntaxRule[] eventrateRules() {
				return new XMLSyntaxRule[]{
				new ElementRule(Parameter.class)
		};
	}


    private XMLSyntaxRule[] popscalingfactorRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class)
        };
    }


    private XMLSyntaxRule[] tippopmodelRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ParametricDistributionModel.class)
        };
    }



    private XMLSyntaxRule[] rootpopmodelRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ParametricDistributionModel.class)
        };
    }


    private XMLSyntaxRule[] hybpopmodelRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ParametricDistributionModel.class)
        };
    }


	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				XMLUnits.SYNTAX_RULES[0],
                new ElementRule(EVENTRATE, eventrateRules()),
                new ElementRule(POPULATION_SCALING_FACTOR, popscalingfactorRules()),
                new ElementRule(TIP_POPULATION_DISTRIBUTION, tippopmodelRules()),
                new ElementRule(ROOT_POPULATION_DISTRIBUTION, rootpopmodelRules()),
                new ElementRule(HYBRID_POPULATION_DISTRIBUTION, hybpopmodelRules())
        };
	}

	@Override
	public String getParserDescription() {
		return "Model for speciation, extinction, hybridization in allopolyploid network.";
	}

	@Override
	public Class getReturnType() {
		return AlloppNetworkPriorModel.class;
	}

}
