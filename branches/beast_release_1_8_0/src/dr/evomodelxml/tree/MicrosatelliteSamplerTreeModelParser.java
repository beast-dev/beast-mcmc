package dr.evomodelxml.tree;

import dr.evolution.alignment.Patterns;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.HashMap;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MicrsatelliteSamplerTreeModel
 *
 */
public class MicrosatelliteSamplerTreeModelParser extends AbstractXMLObjectParser {
    public static final String TREE_MICROSATELLITE_SAMPLER_MODEL = "microsatelliteSamplerTreeModel";
    public static final String TREE = "tree";
    public static final String INTERNAL_VALUES = "internalValues";
    public static final String EXTERNAL_VALUES = "externalValues";
    public static final String USE_PROVIDED_INTERNAL_VALUES = "provideInternalNodeValues";

    public String getParserName(){
        return TREE_MICROSATELLITE_SAMPLER_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel  tree = (TreeModel)xo.getElementFirstChild(TREE);
        Parameter internalVal = (Parameter)xo.getElementFirstChild(INTERNAL_VALUES);
        Patterns microsatellitePattern = (Patterns)xo.getElementFirstChild(EXTERNAL_VALUES);
        int[] externalValues = microsatellitePattern.getPattern(0);
        HashMap<String, Integer> taxaMap = new HashMap<String, Integer>(externalValues.length);
        boolean internalValuesProvided = false;
        if(xo.hasAttribute(USE_PROVIDED_INTERNAL_VALUES)){
            internalValuesProvided = xo.getBooleanAttribute(USE_PROVIDED_INTERNAL_VALUES);
        }
        for(int i = 0; i < externalValues.length; i++){
            taxaMap.put(microsatellitePattern.getTaxonId(i),i);
        }

        String modelId = xo.getAttribute("id", "treeMicrosatelliteSamplerModel");
        return new MicrosatelliteSamplerTreeModel(modelId, tree, internalVal, microsatellitePattern, externalValues, taxaMap, internalValuesProvided);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }
    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_PROVIDED_INTERNAL_VALUES, true),
            new ElementRule(TREE, TreeModel.class),
            new ElementRule(INTERNAL_VALUES, Parameter.class),
            new ElementRule(EXTERNAL_VALUES, Patterns.class)
    };

    public String getParserDescription(){
        return "This parser returns a TreeMicrosatelliteSamplerModel Object";
    }

    public Class getReturnType(){
        return MicrosatelliteSamplerTreeModel.class;

    }
}
