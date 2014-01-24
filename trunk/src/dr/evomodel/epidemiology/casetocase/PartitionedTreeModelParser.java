package dr.evomodel.epidemiology.casetocase;

import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

/**
 * Created by mhall on 24/01/2014.
 */
public class PartitionedTreeModelParser extends TreeModelParser {

    public String getParserName() {
        return PartitionedTreeModel.PARTITIONED_TREE_MODEL;
    }

    // todo Eventually, initialising the branch map should be done here

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel basicTreeModel = (TreeModel)super.parseXMLObject(xo);

        return new PartitionedTreeModel(basicTreeModel);

    }

}