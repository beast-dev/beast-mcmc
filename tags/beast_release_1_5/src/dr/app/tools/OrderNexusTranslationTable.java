package dr.app.tools;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;

import java.io.IOException;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.text.NumberFormat;

/**
 * A simple routine to alphabetically order the taxa translation table
 * in a nexus trees file.
 *
 * @author Alexei Drummond
 */
public class OrderNexusTranslationTable {

    public static void main(String[] args) throws Importer.ImportException, IOException {

        NexusImporter importer = new NexusImporter(new FileReader(args[0]));

        Tree[] trees = importer.importTrees(null);
        System.out.println("Read " + trees.length + " trees from " + args[0]);

        String newFileName = args[0] + ".new";

        PrintStream ps = new PrintStream(new FileOutputStream(newFileName));

        NexusExporter exporter = new NexusExporter(ps);

        exporter.setTreePrefix("STATE_");
        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(7);
        exporter.setNumberFormat(format);
        exporter.setSortedTranslationTable(true);

        exporter.exportTrees(trees, false);

        ps.flush();
        ps.close();
        System.out.println("Wrote " + trees.length + " trees to " + newFileName);
    }

}
