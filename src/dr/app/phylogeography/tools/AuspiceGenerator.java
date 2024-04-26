package dr.app.phylogeography.tools;

import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.SortedRootedTree;
import org.apache.commons.cli.*;

import javax.json.Json;
import javax.json.stream.JsonGenerator;
import javax.json.stream.JsonGeneratorFactory;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDate;
import java.time.Year;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class AuspiceGenerator {
    enum SortType {
        OFF,
        INCREASING,
        DECREASING
    };

    private final double youngestDate;
    private final String[] attributeNames;
    private final Map<String, Set<String>> attributeValues = new HashMap<>();

    public AuspiceGenerator(String treeFileName, String[] attributeNames, SortType sortType, boolean generateMetadata, String outputStem) {

        youngestDate = 2019.0;

        this.attributeNames = attributeNames;

        RootedTree tree;

        try {
            NexusImporter importer = new NexusImporter(new FileReader(treeFileName));
            tree = (RootedTree)importer.importNextTree();
        } catch (ImportException | IOException e) {
            e.printStackTrace();
            return;
        }

        if (sortType == SortType.INCREASING) {
            tree = new SortedRootedTree(tree, SortedRootedTree.BranchOrdering.INCREASING_NODE_DENSITY);
        } else if (sortType == SortType.DECREASING) {
            tree = new SortedRootedTree(tree, SortedRootedTree.BranchOrdering.DECREASING_NODE_DENSITY);
        }

        try {
            writeTreeJSON(outputStem + "_tree.json", tree);

            if (generateMetadata) {
                writeMetaJSON(outputStem + "_meta.json", tree);
            }
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
    }

    private double currentYValue = 0;
    private int cladeNumber = 1;

    private void writeTreeJSON(String filename, RootedTree tree) throws IOException {
        FileWriter writer = new FileWriter(filename);

        HashMap<String, Object> config = new HashMap<>();
        config.put(JsonGenerator.PRETTY_PRINTING, true);

        JsonGeneratorFactory factory = Json.createGeneratorFactory(config);

        JsonGenerator generator = factory.createGenerator(writer);
        currentYValue = 0;

        writeNode(generator, tree, tree.getRootNode());

        generator.close();
    }

    private double writeNode(JsonGenerator generator, RootedTree tree, Node node) {

        double yValue;

        generator.writeStartObject();

        generator
                .writeStartObject("attr")
                .write("clock_length", tree.getLength(node));

        for (String attribute : attributeNames) {
            if (node.getAttribute(attribute) != null) {
                String value = node.getAttribute(attribute).toString();
                generator.write(attribute, value);
                addAttributeValue(attribute, value);

                if (node.getAttribute(attribute + ".set") != null) {
                    Object[] attributeSet = (Object[])node.getAttribute(attribute + ".set");
                    Object[] attributeProbs = (Object[])node.getAttribute(attribute + ".set.prob");

                    generator.writeStartObject(attribute + "_confidence");

                    int i = 0;
                    for (Object attr :  attributeSet) {
                        addAttributeValue(attribute, attr.toString());
                        generator.write(attr.toString(), Double.parseDouble(attributeProbs[i].toString()));
                        i++;
                    }
                    generator.writeEnd();

                } else  {
                    generator
                            .writeStartObject(attribute + "_confidence")
                            .write(node.getAttribute(attribute).toString(), 1.0)
                            .writeEnd();
                }
            } else {
                if (tree.isExternal(node)) {
                    System.out.println("Attribute, " + attribute + ", missing for tip " + tree.getTaxon(node).getName());
                }
            }
        }

        double numDate = youngestDate - tree.getHeight(node);

        double numDateUpper = Double.NaN;
        double numDateLower = Double.NaN;
        if (node.getAttribute("height_95%_HPD") != null) {
            Object[] values = (Object[]) node.getAttribute("height_95%_HPD");
            numDateLower = youngestDate - Double.parseDouble(values[0].toString());
            numDateUpper = youngestDate - Double.parseDouble(values[1].toString());
        }

        String strain;

        String date = numberToCalendarDate(numDate);

        if (!tree.isExternal(node)) {
            strain = "0";
        } else {
//            String[] parts = tree.getTaxon(node).getName().split("\\|");
//            date = parts[parts.length - 1];

            strain =  tree.getTaxon(node).getName();

        }
        generator
                .write("date", date)
//                .write("raw_date", date)
                .write("num_date", numDate);
        if (!Double.isNaN(numDateLower)) {
            generator.writeStartArray("num_date_confidence").write(numDateUpper).write(numDateLower).writeEnd();
        }
        generator
                .write("div", 0);
        generator.writeEnd(); // end "attr"


        if (!tree.isExternal(node)) {
//            generator
//                    .write("date", node.getAttribute("date").toString())
//                    .write("num_date", tree.getHeight(node))
//                    .writeEnd();

            generator.writeStartArray("children");

            yValue = 0.0;
            for (Node child : tree.getChildren(node)) {
                yValue += writeNode(generator, tree, child);
            }
            yValue /= tree.getChildren(node).size();

            generator.writeEnd(); // end "children"

        } else {
            yValue = currentYValue;
            currentYValue += 1.0;
        }

        double branchLength = tree.getLength(node);
        if ( node.getAttribute("rate") != null) {
            branchLength *=  Double.parseDouble(node.getAttribute("rate").toString());
        }
        generator
                .write("branch_length", branchLength)
                .write("strain", strain)
                .write("clade", cladeNumber)
                .write("yvalue", yValue);

        cladeNumber += 1;

        generator.writeEnd(); // object

        return yValue;
    }

    private void addAttributeValue(String attribute, String value) {
        // add the value to the list of possible values
        Set<String> valueSet = attributeValues.computeIfAbsent(attribute, k -> new HashSet<>());
        valueSet.add(value);
    }

    private void writeMetaJSON(String filename, RootedTree tree) throws IOException {
        FileWriter writer = new FileWriter(filename);

        HashMap<String, Object> config = new HashMap<>();
        config.put(JsonGenerator.PRETTY_PRINTING, true);

        JsonGeneratorFactory factory = Json.createGeneratorFactory(config);

        JsonGenerator generator = factory.createGenerator(writer);

        generator.writeStartObject();

        generator.writeStartObject("author_info").writeEnd();

        generator.writeStartObject("color_options");

        for (String attribute : attributeNames) {
            generator.writeStartObject(attribute);

            generator.writeStartArray("color_map");
            for (String value : attributeValues.get(attribute)) {
                generator
                        .writeStartArray()
                        .write(value)
                        .write("#000000")
                        .writeEnd();
            }
            generator.writeEnd(); // color_map array

            generator
                    .write("key", attribute)
                    .write("legendTitle", attribute)
                    .write("menuItem", attribute)
                    .write("type", "discrete");

            generator.writeEnd(); // attribute object
        }

        generator
                .writeStartObject("num_date")
                .write("key", "num_date")
                .write("legendTitle", "Sampling Date")
                .write("menuItem", "date")
                .write("type", "continuous")
                .writeEnd(); // num_date object

        generator.writeEnd(); // color_options object

        generator
                .writeStartObject("defaults")
                .write("mapTriplicate", true)
                .writeEnd(); // defaults object

        generator.writeStartArray("filters");
        for (String attribute : attributeNames) {
            generator.write("attributes");
        }
        generator.writeEnd(); // filters array

        generator.writeStartObject("geo");
        for (String attribute : attributeNames) {
            generator.writeStartObject(attribute);

            for (String value : attributeValues.get(attribute)) {
                generator.writeStartObject(value)
                        .write("latitude", 0.0)
                        .write("longitude", 0.0)
                        .writeEnd();
            }

            generator.writeEnd(); // attribute object
        }
        generator.writeEnd(); // geo object

        generator
                .writeStartArray("maintainer")
                .writeEnd(); // maintainer array

        generator
                .writeStartArray("panels")
                .write("tree")
                .write("map")
                .writeEnd(); // panels array

        generator
                .write("title", "")
                .write("updated", LocalDate.now().toString())
                .write("virus_count", tree.getExternalNodes().size());

        generator.writeEnd(); // outer object

        generator.close();
    }

    private String numberToCalendarDate(double numDate) {
        return Year.of( (int)Math.floor(numDate) )
                .atDay( (int)Math.floor((numDate - Math.floor(numDate)) * 365) + 1 )
                .format(DateTimeFormatter.ISO_LOCAL_DATE);
    }

    public static void main(String[] args) {
        // create Options object
        Options options = new Options();

        options.addOption("h", "help", false, "display help");
        options.addOption("v", "verbose", false, "write analysis details to stderr");

        options.addOption("m", "metadata", false, "generate a metadata file");

        options.addOption( Option.builder( "s" )
                .longOpt("sort")
                .argName("up/down")
                .hasArg()
                .required(false)
                .desc( "sort the nodes in increasing node density (up) or decreasing (down)" )
                .type(String.class).build());

        options.addOption( Option.builder( "t" )
                .longOpt("tree")
                .argName("file name")
                .hasArg()
                .required(true)
                .desc( "input tree file" )
                .type(String.class).build());

        options.addOption( Option.builder( "a" )
                .longOpt("attributes")
                .argName("name")
                .hasArg()
                .numberOfArgs(Option.UNLIMITED_VALUES)
                .required(true)
                .desc( "a list of attribute names" )
                .type(String.class).build());


        options.addOption( Option.builder( "o" )
                .longOpt("output")
                .argName("file name")
                .hasArg()
                .required(true)
                .desc( "output tree file" )
                .type(String.class).build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse( options, args);
        } catch (ParseException pe) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "beast2auspice", options, true );
            return;
        }

        if(cmd.hasOption("h")) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "beast2auspice", options, true );
            return;
        }

        boolean verbose = cmd.hasOption("v");

        boolean generateMetadata = cmd.hasOption("m");

        SortType sortType = SortType.OFF;

        if (cmd.hasOption("s")) {
            String sort = cmd.getOptionValue("s");
            if (sort.equalsIgnoreCase("up")) {
                sortType = SortType.INCREASING;
            } else if (sort.equalsIgnoreCase("down")) {
                sortType = SortType.DECREASING;
            } else {
                System.err.println("Unrecognised sort option: " + sort);
            }
        }

        String treeFileName = cmd.getOptionValue("t");
        String outputFilename = cmd.getOptionValue("o");

        String[] attributeNames = cmd.getOptionValues("a");

        if (verbose) {
            System.err.println("     input tree file: " + treeFileName);
            System.err.println("output filename stem: " + outputFilename);
        }

        new AuspiceGenerator(treeFileName, attributeNames, sortType, generateMetadata, outputFilename);

    }

}
