package dr.evomodelxml.speciation;

import dr.evomodel.speciation.BirthDeathCollapseModel;
import dr.evomodel.speciation.BirthDeathCollapseNClustersStatistic;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.xml.*;

/**
 * @author Graham Jones
 *         Date: 01/09/2013
 */
public class BirthDeathCollapseNClustersStatisticParser extends AbstractXMLObjectParser {
    public static final String BDC_NCLUSTERS_STATISTIC = "bdcNClustersStatistic";
    public static final String SPECIES_TREE = "speciesTree";
    public static final String COLLAPSE_MODEL = "collapseModel";



    public String getParserName() {
        return BDC_NCLUSTERS_STATISTIC;
    }


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        xo.getAttribute("name");
        final XMLObject spptreexo = xo.getChild(SPECIES_TREE);
        SpeciesTreeModel spptree = (SpeciesTreeModel) spptreexo.getChild(SpeciesTreeModel.class);
        final XMLObject cmxo = xo.getChild(COLLAPSE_MODEL);
        BirthDeathCollapseModel bdcm = (BirthDeathCollapseModel)cmxo.getChild(BirthDeathCollapseModel.class);
        return new BirthDeathCollapseNClustersStatistic(spptree, bdcm);
    }


    private  XMLSyntaxRule[] spptreeRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciesTreeModel.class)
        };
    }

    private  XMLSyntaxRule[] bdcmRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(BirthDeathCollapseModel.class)
        };
    }



    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newStringRule("name"),
                new ElementRule(SPECIES_TREE, spptreeRules()),
                new ElementRule(COLLAPSE_MODEL, bdcmRules())
        };
    }


    @Override
    public String getParserDescription() {
        return "Statistic for number of collapsed nodes in species tree when using birth-death-collapse model.";
    }

    @Override
    public Class getReturnType() {
        return BirthDeathCollapseNClustersStatistic.class;
    }


}
