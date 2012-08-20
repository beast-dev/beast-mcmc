package dr.evomodelxml.speciation;


import dr.util.Attributable;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;


/**
 * 
 * Parses a AlloppSpeciesNetworkModel.
 *
 * @author Graham Jones
 *         Date: 19/04/2011
 */



/*
 * 
 * Parses a MulLabSpeciesTreeModel which might look like this.
 *
  <alloppSpeciesNetwork id="apspnetwork" oneHybridization="true">
    <alloppSpecies idref="alloppSpecies"/>
    <sppSplitPopulations value="0.018">
      <parameter id="apspnetwork.splitPopSize"/>
    </sppSplitPopulations>
    <preRootHeights>
      <parameter id="pre.root.heights"/>
    </preRootHeights>
  </alloppSpeciesNetwork>
 *
 */


/* 
 * 2011-07-07 decided to remove preRootHeights, and therefore these bits of code.
 * (instead, incorporate into AlloppNetworkNodeSlide)
 * 
 * 
 * public static final String PRE_ROOT_HEIGHTS = "preRootHeights";
 * 
 * 
 * final XMLObject prhxo = xo.getChild(PRE_ROOT_HEIGHTS);
 * 
 * 
 *  ParameterParser.replaceParameter(prhxo, asnm.prHeights);
 *  final Parameter.DefaultBounds prhbounds =
 *  new Parameter.DefaultBounds(Double.MAX_VALUE, 0, asnm.prHeights.getDimension());
 *  asnm.prHeights.addBounds(prhbounds);
 * 
 * 	private ElementRule prhElementRule() {
 * 	return new ElementRule(PRE_ROOT_HEIGHTS, new XMLSyntaxRule[]{
 *                         new ElementRule(Parameter.class)});
 *                         }
 *                         
 *    also code in    AlloppSpeciesNetworkModel                  
 *    
 */



public class AlloppSpeciesNetworkModelParser extends AbstractXMLObjectParser {
	public static final String ALLOPPSPECIESNETWORK = "alloppSpeciesNetwork";
    public static final String ONEHYBRIDIZATION = "oneHybridization";
    public static final String SPP_SPLIT_POPULATIONS = "sppSplitPopulations";


	
	
	public String getParserName() {
		return ALLOPPSPECIESNETWORK;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		AlloppSpeciesBindings apspb = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
		boolean onehyb = xo.getBooleanAttribute(ONEHYBRIDIZATION);
		
		final XMLObject sppxo = xo.getChild(SPP_SPLIT_POPULATIONS);
		final double sppvalue = sppxo.getAttribute(Attributable.VALUE, 1.0);
		
		// JH in SpeciesTreeModelParser creates Parameter sppSplitPopulations
		// here, calling SpeciesTreeModel.createSplitPopulationsParameter()
		// which calls SpeciesBindings.nSpecies(). Seems complicated.
        // Something must be done (replaceParameter) or MixedDistributionLikelihood
		// gives an error data length != indicator length.
		AlloppSpeciesNetworkModel asnm = new AlloppSpeciesNetworkModel(apspb, sppvalue, onehyb);
		
        ParameterParser.replaceParameter(sppxo, asnm.popvalues);
        final Parameter.DefaultBounds sppbounds =
                new Parameter.DefaultBounds(Double.MAX_VALUE, 0, asnm.popvalues.getDimension());
        asnm.popvalues.addBounds(sppbounds);
        
		return asnm;
	}

	
	private ElementRule sppElementRule() {
		return new ElementRule(SPP_SPLIT_POPULATIONS, new XMLSyntaxRule[]{
	                        AttributeRule.newDoubleRule(Attributable.VALUE, true),
	                        new ElementRule(Parameter.class)});
		}
	
	

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
        		AttributeRule.newBooleanRule(ONEHYBRIDIZATION, true),
                new ElementRule(AlloppSpeciesBindings.class),
                sppElementRule()
        };
	}

	@Override
	public String getParserDescription() {
		return "Species network with population sizes along branches";
	}

	@Override
	public Class getReturnType() {
		return AlloppSpeciesNetworkModel.class;
	}

}
