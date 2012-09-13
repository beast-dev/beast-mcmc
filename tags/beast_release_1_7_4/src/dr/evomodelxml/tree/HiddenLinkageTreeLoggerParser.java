package dr.evomodelxml.tree;

import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.HiddenLinkageTreeLogger;
import dr.evomodel.tree.TreeLogger;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling (koadman)
 */
public class HiddenLinkageTreeLoggerParser extends TreeLoggerParser {

    private XMLSyntaxRule[] rules;

	public HiddenLinkageTreeLoggerParser(){
		rules = new XMLSyntaxRule[super.getSyntaxRules().length + 1];
		System.arraycopy(super.getSyntaxRules(), 0, rules, 0, rules.length-1);
		rules[rules.length-1] = new ElementRule(HiddenLinkageModel.class);
	}
    public String getParserName() {
        return "logHiddenLinkageTree";
    }

    public String getParserDescription() {
        return "Logs a tree with hidden linkage among metagenomic reads to a file";
    }

    public Class getReturnType() {
        return HiddenLinkageTreeLogger.class;
    }
    
    public XMLSyntaxRule[] getSyntaxRules() {
        return this.rules;
    }
    
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	parseXMLParameters(xo);
    	HiddenLinkageModel hlm = (HiddenLinkageModel)xo.getChild(HiddenLinkageModel.class);

    	HiddenLinkageTreeLogger logger = new HiddenLinkageTreeLogger(hlm, tree, branchRates,
                treeAttributeProviders, treeTraitProviders,
                formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, condition/*,
                normaliseMeanRateTo*/);

        if (title != null) {
            logger.setTitle(title);
        }

        return logger;
    }
}
