package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 */
public class NodeTraitLogger extends TreeTraitProvider.Helper {

	public static final String TRAIT_LOGGER = "logAllTraits";

	private final TreeModel treeModel;

	public NodeTraitLogger(TreeModel treeModel) {
		this.treeModel = treeModel;

        throw new RuntimeException("NodeTraitLogger has not been fully converted to using TreeTraitProvider");

        // need to find all the traits here....
	}

	public static String[] getAllNodeTraitLabels(TreeModel tree) {

		Map<String, Parameter> traits = tree.getTraitMap(tree.getRoot());
		List<String> labels = new ArrayList<String>();
		for (Map.Entry<String, Parameter> stringParameterEntry : traits.entrySet()) {

			Parameter traitParameter = stringParameterEntry.getValue();
			if (traitParameter.getDimension() == 1)
				labels.add(stringParameterEntry.getKey());
			else {
				for (int i = 1; i <= traitParameter.getDimension(); i++)
					labels.add(stringParameterEntry.getKey() + i);
			}
		}
		return labels.toArray(new String[labels.size()]);
	}

	public static String[] getAllNodeTraitValues(TreeModel tree, NodeRef node) {

		Map<String, Parameter> traits = tree.getTraitMap(node);
		List<String> values = new ArrayList<String>();
		for( Map.Entry<String, Parameter> stringParameterEntry : traits.entrySet() ) {

			Parameter traitParameter = stringParameterEntry.getValue();

			for (int i = 0; i < traitParameter.getDimension(); i++)
				values.add(Double.toString(traitParameter.getParameterValue(i)));
		}
		return values.toArray(new String[values.size()]);
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return TRAIT_LOGGER;
		}

		/**
		 * @return an object based on the XML element it was passed.
		 */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

			try {

				//Map<String, Parameter> traits = treeModel.getTraitMap(treeModel.getRoot());
				return new NodeTraitLogger(treeModel);

			} catch (IllegalArgumentException e) {

				throw new XMLParseException("Tree " + treeModel.getId() + " contains no traits to log");
			}


		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				new ElementRule(TreeModel.class, "The tree which is to be logged")
		};

		public String getParserDescription() {
			return null;
		}

		public String getExample() {
			return null;
		}

		public Class getReturnType() {
			return TreeTraitProvider.class;
		}
	};
}
