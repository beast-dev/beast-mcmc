package dr.app.beauti.components.dnds;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.util.XMLWriter;
import dr.util.Attribute;

/**
 * @author Filip Bielejec
 * @version $Id$
 */

public class DnDsComponentGenerator extends BaseComponentGenerator {

	protected DnDsComponentGenerator(BeautiOptions options) {
		super(options);
	}

	public boolean usesInsertionPoint(InsertionPoint point) {

		 DnDsComponentOptions component = (DnDsComponentOptions)
                 options.getComponentOptions(DnDsComponentOptions.class);

        if (component.getPartitionList().size() == 0) {
            // Empty, so do nothing
            return false;
        }

		switch (point) {
		case IN_OPERATORS:
		case IN_FILE_LOG_PARAMETERS:
		case IN_TREES_LOG:
		case AFTER_TREES_LOG:
			return true;
		default:
			return false;
		}

	}

	protected void generate(InsertionPoint point, Object item, XMLWriter writer) {

		DnDsComponentOptions component = (DnDsComponentOptions) options
				.getComponentOptions(DnDsComponentOptions.class);

		switch (point) {
		case IN_OPERATORS:
			writeCodonPartitionedRobustCounting(writer, component);
			break;
		case IN_FILE_LOG_PARAMETERS:
			// TODO
			break;
		case IN_TREES_LOG:
			// TODO
			break;
		case AFTER_TREES_LOG:
			break;
		default:
			throw new IllegalArgumentException(
					"This insertion point is not implemented for "
							+ this.getClass().getName());
		}

	}

	protected String getCommentLabel() {
		return "Codon partitioned robust counting";
	}

	private void writeCodonPartitionedRobustCounting(XMLWriter writer,
			DnDsComponentOptions component) {

        for (PartitionSubstitutionModel model : component.getPartitionList()) {
            writeCodonPartitionedRobustCounting(writer, model);
        }
	}// END: writeCodonPartitionedRobustCounting()

    // Called for each model that requires robust counting (can be more than one)
    private void writeCodonPartitionedRobustCounting(XMLWriter writer,
                                                     PartitionSubstitutionModel model) {

        System.err.println("DEBUG: Writing RB for " + model.getName());
        writer.writeComment("Robust counting for: " + model.getName());
//		// Just sth to see if it works at all
//		writer.writeOpenTag("codonPartitionedRobustCounting",
//				new Attribute[] {
//						new Attribute.Default<String>("id", "robustCounting1"),
//						new Attribute.Default<String>("labeling", "S"),
//						new Attribute.Default<String>("useUniformization",
//								"true"),
//						new Attribute.Default<String>("unconditionedPerBranch",
//								"true") });
//
//		writer.writeIDref("treeModel", "");
//		writer.writeCloseTag("codonPartitionedRobustCounting");
    }

}
