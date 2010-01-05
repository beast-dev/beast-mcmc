package dr.app.phylogeography;

import dr.app.phylogeography.spread.Location;
import dr.app.phylogeography.structure.Layer;
import dr.app.phylogeography.structure.TimeLine;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;

import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class Phylogeography {
    // used to calibrate the color range for the branches
    double minRate;
    double maxRate;
    double medianRate;

    double maxHeight;
    double maxBranchLength;

    private final RootedTree tree;
    private final Map<String, Location> locationMap;

    private List<Layer> layers = new ArrayList<Layer>();

    public Phylogeography(RootedTree tree) {
        this(tree, null);
    }

    public Phylogeography(RootedTree tree, List<Location> locationList) {

        this.tree = tree;
        maxHeight = tree.getHeight(tree.getRootNode());
        maxBranchLength = 0.0;

        for (Node node : tree.getNodes()) {
            if (!tree.isRoot(node)) {
                if (tree.getLength(node) > maxBranchLength) {
                    maxBranchLength = tree.getLength(node);
                }
            }
        }

        if (locationList != null) {
            double minLat = Double.MAX_VALUE;
            double maxLat = -Double.MAX_VALUE;
            double minLong = Double.MAX_VALUE;
            double maxLong = -Double.MAX_VALUE;

            locationMap = new HashMap<String, Location>();

            for (Location location : locationList) {
                if (location.getLatitude() < minLat)  {
                    minLat = location.getLatitude();
                }
                if (location.getLatitude() > maxLat)  {
                    maxLat = location.getLatitude();
                }

                if (location.getLongitude() < minLong)  {
                    minLong = location.getLongitude();
                }
                if (location.getLongitude() > maxLong)  {
                    maxLong = location.getLongitude();
                }
                locationMap.put(location.getState(), location);
            }

//            radius = 100*Math.abs(maxLat-minLat)*Math.abs(maxLong-minLong);
//            radius = 200000;

//            if (mostRecentDate - treeToExport.getHeight(treeToExport.getRootNode()) < 0) {
//                ancient = true;
//            }

        } else {
            locationMap = null;
        }

    }

    public List<Layer> getLayers() {
        return layers;
    }

    public static void main(String[] args) {

        String inputTreeFile = args[0];
        String outputFile = args[1];
        RootedTree tree = null;

        try {
            TreeImporter importer = new NexusImporter(new FileReader(inputTreeFile));
            tree = (RootedTree)importer.importNextTree();
        } catch (
                ImportException e) {
            e.printStackTrace();
            return;
        } catch (
                IOException e) {
            e.printStackTrace();
            return;
        }

        Phylogeography phylogeography = new Phylogeography(tree);

        Collection<Layer> layers = phylogeography.getLayers();

        TimeLine timeLine = new TimeLine(2000, 2009, 100);
        
//        try {
//            Generator generator = new KMLGenerator(timeLine);
//
//            generator.generate(new PrintWriter(outputFile), layers);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }

    }
}
