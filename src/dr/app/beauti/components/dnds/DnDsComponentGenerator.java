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

		DnDsComponentOptions component = (DnDsComponentOptions) options
				.getComponentOptions(DnDsComponentOptions.class);

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

	}// END: usesInsertionPoint

	protected void generate(InsertionPoint point, Object item, XMLWriter writer) {

		DnDsComponentOptions component = (DnDsComponentOptions) options
				.getComponentOptions(DnDsComponentOptions.class);

		// writePatternList(writer, component);

		switch (point) {
		case IN_OPERATORS:
			writeCodonPartitionedRobustCounting(writer, component);
			break;
		case IN_FILE_LOG_PARAMETERS:
			writeLogs(writer, component);
			break;
		case IN_TREES_LOG:
			writeTreeLogs(writer, component);
			break;
		case AFTER_TREES_LOG:
			writeDNdSLogger(writer, component);
			writeDNdSReport(writer, component);
			break;
		default:
			throw new IllegalArgumentException(
					"This insertion point is not implemented for "
							+ this.getClass().getName());
		}

	}// END: generate

	protected String getCommentLabel() {
		return "Codon partitioned robust counting";
	}

	private void writeCodonPartitionedRobustCounting(XMLWriter writer,
			DnDsComponentOptions component) {

		for (PartitionSubstitutionModel model : component.getPartitionList()) {
			writeCodonPartitionedRobustCounting(writer, model);
		}
	}

	// Called for each model that requires robust counting (can be more than
	// one)
	private void writeCodonPartitionedRobustCounting(XMLWriter writer,
			PartitionSubstitutionModel model) {

		System.err.println("DEBUG: Writing RB for " + model.getName());

		writer.writeComment("Robust counting for: " + model.getName());

		// TODO: Hand coding is so 90s

		// S operator
		writer.writeOpenTag("codonPartitionedRobustCounting",
				new Attribute[] {
						new Attribute.Default<String>("id", "robustCounting1"),
						new Attribute.Default<String>("labeling", "S"),
						new Attribute.Default<String>("useUniformization",
								"true"),
						new Attribute.Default<String>("unconditionedPerBranch",
								"true") });

		writer.writeIDref("idref", "treeModel");
		writer.writeOpenTag("firstPosition");
		writer.writeIDref("ancestralTreeLikelihood", "CP1.treeLikelihood");
		writer.writeCloseTag("firstPosition");

		writer.writeOpenTag("secondPosition");
		writer.writeIDref("ancestralTreeLikelihood", "CP2.treeLikelihood");
		writer.writeCloseTag("secondPosition");

		writer.writeOpenTag("thirdPosition");
		writer.writeIDref("ancestralTreeLikelihood", "CP3.treeLikelihood");
		writer.writeCloseTag("thirdPosition");

		writer.writeCloseTag("codonPartitionedRobustCounting");

		writer.writeBlankLine();

		// N operator:
		writer.writeOpenTag("codonPartitionedRobustCounting",
				new Attribute[] {
						new Attribute.Default<String>("id", "robustCounting2"),
						new Attribute.Default<String>("labeling", "N"),
						new Attribute.Default<String>("useUniformization",
								"true"),
						new Attribute.Default<String>("unconditionedPerBranch",
								"true") });

		writer.writeIDref("idref", "treeModel");
		writer.writeOpenTag("firstPosition");
		writer.writeIDref("ancestralTreeLikelihood", "CP1.treeLikelihood");
		writer.writeCloseTag("firstPosition");

		writer.writeOpenTag("secondPosition");
		writer.writeIDref("ancestralTreeLikelihood", "CP2.treeLikelihood");
		writer.writeCloseTag("secondPosition");

		writer.writeOpenTag("thirdPosition");
		writer.writeIDref("ancestralTreeLikelihood", "CP3.treeLikelihood");
		writer.writeCloseTag("thirdPosition");

		writer.writeCloseTag("codonPartitionedRobustCounting");

	}// END: writeCodonPartitionedRobustCounting()

	private void writeLogs(XMLWriter writer, DnDsComponentOptions component) {

		for (PartitionSubstitutionModel model : component.getPartitionList()) {
			writeLogs(writer, model);
		}
	}

	private void writeLogs(XMLWriter writer, PartitionSubstitutionModel model) {

		writer.writeComment("Robust counting for: " + model.getName());

		writer.writeIDref("codonPartitionedRobustCounting", "robustCounting1");
		writer.writeIDref("codonPartitionedRobustCounting", "robustCounting2");

	}// END: writeLogs()

	private void writeTreeLogs(XMLWriter writer, DnDsComponentOptions component) {

		for (PartitionSubstitutionModel model : component.getPartitionList()) {
			// We re-use the same method
			writeLogs(writer, model);
		}
	}

	private void writeDNdSLogger(XMLWriter writer,
			DnDsComponentOptions component) {

		for (PartitionSubstitutionModel model : component.getPartitionList()) {
			writeDNdSLogger(writer, model);
		}
	}

	private void writeDNdSLogger(XMLWriter writer,
			PartitionSubstitutionModel model) {

		writer.writeComment("Robust counting for: " + model.getName());

		writer.writeOpenTag("log", new Attribute[] {
				new Attribute.Default<String>("id", "fileLog_dNdS"),
				new Attribute.Default<String>("logEvery", "10000"),
				new Attribute.Default<String>("fileName", model.getName()
						+ ".log") });

		writer
				.writeOpenTag("dNdSLogger",
						new Attribute[] { new Attribute.Default<String>("id",
								"dNdS") });
		writer.writeIDref("idref", "treeModel");
		writer.writeIDref("codonPartitionedRobustCounting", "robustCounting1");
		writer.writeIDref("codonPartitionedRobustCounting", "robustCounting2");
		writer.writeCloseTag("dNdSLogger");

		writer.writeCloseTag("log");

	}// END: writeLogs()

	private void writeDNdSReport(XMLWriter writer,
			DnDsComponentOptions component) {

		for (PartitionSubstitutionModel model : component.getPartitionList()) {
			writeDNdSReport(writer, model);
		}
	}

	private void writeDNdSReport(XMLWriter writer,
			PartitionSubstitutionModel model) {

		writer.writeComment("Robust counting for: " + model.getName());

		writer.writeOpenTag("report");

		writer.write("<dNdSPerSiteAnalysis fileName=" + '\"' + model.getName()
				+ ".log\"/> \n");

		writer.writeCloseTag("report");

	}// END: writeDNdSReport()

	// private void writePatternList(XMLWriter writer,
	// DnDsComponentOptions component) {
	//
	// try {
	//
	// for (PartitionSubstitutionModel model : component
	// .getPartitionList()) {
	//
	// writer
	// .writeComment("unique=\"false\" returns one pattern for each and every site");
	//
	// // TODO do we need this loop?
	// for (AbstractPartitionData partition : options.dataPartitions) {
	//
	// for (int i = 1; i <= 3; i++) {
	//
	// writePatternList((PartitionData) partition, i, 3,
	// false, model.getPrefix(i), writer);
	//
	// writer.writeBlankLine();
	//
	// }// END: codon positions loop
	// }// END: AbstractPartitionData loop
	//
	// }// END: model loop
	//
	// } catch (Exception e) {
	// e.getCause().printStackTrace();
	// }
	//
	// }// END: modifySitePatterns()

	// @SuppressWarnings("unchecked")
	// private void writePatternList(PartitionData partition, int from, int
	// every,
	// boolean unique, String codonPrefix, XMLWriter writer) {
	//
	// Alignment alignment = partition.getAlignment();
	//
	// List<Attribute> attributes = new ArrayList<Attribute>();
	//
	// attributes.add(new Attribute.Default<String>(XMLParser.ID, codonPrefix
	// + partition.getPrefix() + SitePatternsParser.PATTERNS));
	//
	// attributes.add(new Attribute.Default<String>("from", "" + from));
	//
	// attributes.add(new Attribute.Default<String>("every", "" + every));
	//
	// attributes.add(new Attribute.Default<String>("unique", "" + unique));
	//
	// writer.writeOpenTag(SitePatternsParser.PATTERNS, attributes);
	//
	// writer.writeIDref(AlignmentParser.ALIGNMENT, alignment.getId());
	//
	// writer.writeCloseTag(SitePatternsParser.PATTERNS);
	//
	// }

}// END: class
