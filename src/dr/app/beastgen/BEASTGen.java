package dr.app.beastgen;

import dr.evolution.io.Importer;
import freemarker.template.*;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class BEASTGen {

    private static final String ARGUMENT_CHARACTER = "-";

    public BEASTGen(Map argumentMap, String templateFileName, String inputFileName, String outputFileName) throws IOException {

        Configuration cfg = new Configuration();
        cfg.setObjectWrapper(new DefaultObjectWrapper());
        cfg.setDefaultEncoding("UTF-8");
        cfg.setTemplateExceptionHandler(TemplateExceptionHandler.DEBUG_HANDLER);
        cfg.setIncompatibleImprovements(new Version(2, 3, 20));  // FreeMarker 2.3.20

        Map root = null;
        try {
            root = constructDataModel(inputFileName);
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

    private Map constructDataModel(String inputFileName) throws IOException, Importer.ImportException {
        DataModelImporter importer = new DataModelImporter();

        Map root = importer.importFromFile(new File(inputFileName));

        return root;
    }

    public static String[] parseArguments(String[] arguments, Map dataModel) {

        int[] optionIndex = new int[arguments.length];
        for (int i = 0; i < optionIndex.length; i++) {
            optionIndex[i] = -1;
        }

        int i = 0;
        for (; i < arguments.length; i++) {
            if (arguments[i].startsWith(ARGUMENT_CHARACTER)) {

                // the first value may be appended to the option label (e.g., '-t1.0'):
                String tag = arguments[i].substring(1);
                if (!tag.isEmpty() && Character.isLetter(tag.charAt(0))) {
                    i++;
                    if (i == arguments.length) {
                        System.err.println("Value for argument, " + tag + ", is missing");
                        System.exit(1);
                    }
                    dataModel.put(tag, arguments[i]);
                } else {
                    System.err.println("Illegal argument, " + tag);
                    System.exit(1);
                }
            } else {
                break;
            }
        }

        String[] leftoverArguments = new String[arguments.length - i];
        for (int j = 0; j < leftoverArguments.length; j++) {
            leftoverArguments[j] = arguments[i + j];
        }

        return leftoverArguments;
    }

    public static void main(String[] args) {

        Map argumentMap = new HashMap();

        String[] arguments = parseArguments(args, argumentMap);
        if (arguments.length < 2 || arguments.length > 3) {
            System.out.println("Usage: beastgen [-tag=value ...] <template_filename> <nexus_filename> [<output_filename>]");
        }

        try {
            new BEASTGen(argumentMap, arguments[0], arguments[1], (arguments.length == 3 ? arguments[2] : null));
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
}
