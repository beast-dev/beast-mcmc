package dr.app.beastgen;

import dr.app.beauti.options.DateGuesser;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import freemarker.template.*;

import java.io.*;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

public class BEASTGen {

    public BEASTGen(DateGuesser guesser, Map argumentMap, String templateFileName, String inputFileName, String outputFileName) throws IOException {

        Configuration cfg = new Configuration();
        cfg.setObjectWrapper(new DefaultObjectWrapper());
        cfg.setDefaultEncoding("UTF-8");
        cfg.setTemplateExceptionHandler(TemplateExceptionHandler.DEBUG_HANDLER);
        cfg.setIncompatibleImprovements(new Version(2, 3, 20));  // FreeMarker 2.3.20

        Map root = null;
        try {
            root = constructDataModel(inputFileName, guesser);
        } catch (Importer.ImportException ie) {
            System.err.println("Error importing file: " + ie.getMessage());
        }

        root.putAll(argumentMap);

        Template template = cfg.getTemplate(templateFileName);

        try {
            Writer out = (outputFileName != null ?
                    new FileWriter(new File(outputFileName)) :
                    new OutputStreamWriter(System.out));

            template.process(root, out);
        } catch (TemplateException te) {
            System.err.println("Error processing template, " + templateFileName + ": " + te.getMessage());
            System.exit(1);
        }
    }

    private Map constructDataModel(String inputFileName, DateGuesser guesser) throws IOException, Importer.ImportException {
        DataModelImporter importer = new DataModelImporter(guesser);

        Map root = importer.importFromFile(new File(inputFileName));

        return root;
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("BEASTGen v1.0, 2013", 60);
        centreLine("BEAST input file generator", 60);
        centreLine("Andrew Rambaut, University of Edinburgh", 60);
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beastgen", "<template-file-name> <input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: beastgen template.beast test.nex test.xml");
        System.out.println("  Example: beastgen -help");
        System.out.println();
    }


    public static void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("date_order", "The order of the date field (negative numbers from last)"),
                        new Arguments.StringOption("date_prefix", "prefix", "A string that is the prefix to the date field"),
                        new Arguments.StringOption("date_regex", "regex", "A string that gives the regular expression to match the date"),
                        new Arguments.StringOption("date_format", "format", "A string that gives the date format for parsing"),
                        new Arguments.Option("date_precision", "Specifies the date is a variable precision yyyy-MM-dd format"),
                        new Arguments.StringOption("D", "\"key=value,key=value...\"", "Properties for exchange in templates"),
                        new Arguments.Option("version", "Print the version and credits and stop"),
                        new Arguments.Option("help", "Print this information and stop"),
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            printTitle();
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printUsage(arguments);
            System.exit(1);
        }

        DateGuesser guesser = new DateGuesser();
//        guesser.guessDates = false;
//        guesser.guessType = DateGuesser.GuessType.ORDER;
//        guesser.fromLast = false;
//        guesser.order = 0;
//        guesser.prefix = null;
//        guesser.regex = null;
//        guesser.offset = 0.0;
//        guesser.unlessLessThan = 0.0;
//        guesser.offset2 = 0.0;
//        guesser.parseCalendarDates = false;
//        guesser.parseCalendarDatesAndPrecision = false;
//        guesser.calendarDateFormat = "yyyy-MM-dd";

        if (arguments.hasOption("date_order")) {
            guesser.guessDates = true;
            int order = arguments.getIntegerOption("date_order");
            if (order < 0) {
                guesser.order = 1 + order;
                guesser.fromLast = true;
            } else if (order > 0) {
                guesser.order = order - 1;
                guesser.fromLast = false;
            } else {
                guesser.order = 0;
                guesser.fromLast = false;
            }
        }

        if (arguments.hasOption("date_prefix")) {
            guesser.guessDates = true;
            guesser.prefix = arguments.getStringOption("date_prefix");
        }

        if (arguments.hasOption("date_regex")) {
            guesser.guessDates = true;
            guesser.regex = arguments.getStringOption("date_regex");
        }

        if (arguments.hasOption("date_format")) {
            guesser.guessDates = true;
            guesser.calendarDateFormat = arguments.getStringOption("date_format");
            guesser.parseCalendarDates = true;
        }

        if (arguments.hasOption("date_precision")) {
            guesser.guessDates = true;
            guesser.parseCalendarDatesAndPrecision = true;
        }

        Map argumentMap = new HashMap();

        if (arguments.hasOption("D")) {
            String properties = arguments.getStringOption("D");
            for (String property : properties.split(",\\s")) {
                String[] keyValue = property.split("=");
                if (keyValue.length != 2) {
                    System.err.println("Properties should take the form: key=value");
                    System.exit(1);
                }
                String key = keyValue[0].trim();
                String value = keyValue[1].trim();
                if (key.isEmpty()) {
                    System.err.println("Properties should take the form: key=value");
                    System.exit(1);
                }
                argumentMap.put(key, value);
            }
        }

        if (arguments.hasOption("help")) {
            printTitle();
            printUsage(arguments);
            System.exit(0);
        }

        if (arguments.hasOption("version")) {
            printTitle();
        }


        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length < 2 || args2.length > 3) {
            System.err.println("Unknown option: " + args2[0]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        try {
            new BEASTGen(guesser, argumentMap, args2[0], args2[1], (args2.length == 3 ? args2[2] : null));
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
}
