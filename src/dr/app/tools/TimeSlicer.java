package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.continuous.SphericalPolarCoordinates;
import dr.util.Version;
import dr.util.HeapSort;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.geo.KernelDensityEstimator2D;
import dr.geo.KMLCoordinates;
import dr.geo.Polygon2D;
import dr.geo.contouring.*;
import dr.inference.trace.TraceDistribution;
import dr.evomodel.tree.TreeModel;
import dr.stats.DiscreteStatistics;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.awt.geom.Point2D;

import org.jdom.Element;
import org.jdom.output.XMLOutputter;
import org.jdom.output.Format;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class TimeSlicer {

    public static final String sep = "\t";
    public static final String PRECISION_STRING = "precision";
    public static final String RATE_STRING = "rate"; // TODO the rate used is supposed to be the relaxed diffusion rate, but relaxed clock models will also result in a "rate attribute", we somehow need to distinguis between them!

    public static final String SLICE_ELEMENT = "slice";
    public static final String REGIONS_ELEMENT = "hpdRegion";
    public static final String TRAIT = "trait";
    public static final String LOCATIONTRAIT = "location";
    public static final String NAME = "name";
    public static final String DENSITY_VALUE = "density";
    public static final String SLICE_VALUE = "time";
    
    public static final String STYLE = "Style";
    public static final String ID = "id";
    public static final String WIDTH = "0.5";
    public static final String startHPDColor = "00F1D6"; //blue=B36600
    public static final String endHPDColor = "00FF00"; //red=0000FF
    public static final String opacity = "9f";

    public static final String BURNIN = "burnin";
    public static final String SLICE = "slice";
    public static final String FILE_SLICE = "fileSlice";
    public static final String MRSD = "mrsd";
    public static final String HELP = "help";
    public static final String NOISE = "noise";
    public static final String IMPUTE = "impute";
    public static final String SUMMARY = "summary";
    public static final String FORMAT = "format";
    public static final String MODE = "mode";
    public static final String NORMALIZATION = "normalization";
    public static final String HPD = "hpd";
    public static final String SDR = "sdr";

    public static final String[] falseTrue = new String[] {"false","true"};


    public TimeSlicer(String treeFileName, int burnin, String[] traits, double[] slices, boolean impute,
                      boolean trueNoise, double mrsd, ContourMode contourMode,
                      Normalization normalize, boolean getSRD) {

        this.traits = traits;
        traitCount = traits.length;
        this.slices = slices;
        sliceCount = 1;
        doSlices = false;
        mostRecentSamplingDate = mrsd;
        this.contourMode = contourMode;
        sdr = getSRD;

        if (slices != null) {
            sliceCount = slices.length;
            doSlices = true;
        }

        values = new ArrayList<List<List<Trait>>>(sliceCount);
        for (int i = 0; i < sliceCount; i++) {
            List<List<Trait>> thisSlice = new ArrayList<List<Trait>>(traitCount);
            values.add(thisSlice);
            for (int j = 0; j < traitCount; j++) {
                List<Trait> thisTraitSlice = new ArrayList<Trait>();
                thisSlice.add(thisTraitSlice);
            }
        }

        try {
            readAndAnalyzeTrees(treeFileName, burnin, traits, slices, impute, trueNoise, normalize);
        } catch (IOException e) {
            System.err.println("Error reading file: " + treeFileName);
            System.exit(-1);
        } catch (Importer.ImportException e) {
            System.err.println("Error parsing trees in file: " + treeFileName);
            System.exit(-1);
        }
        progressStream.println(treesRead+" trees read.");
        progressStream.println(treesAnalyzed+" trees analyzed.");

    }

    public void output(String outFileName, boolean summaryOnly) {
        output(outFileName,summaryOnly,OutputFormat.XML, 0.80, null);
    }

    private Element rootElement;
    private Element documentElement;
    private Element folderElement;

    public void output(String outFileName, boolean summaryOnly, OutputFormat outputFormat, double hpdValue, String sdrFile) {

        resultsStream = System.out;

        if (outFileName != null) {
            try {
            resultsStream = new PrintStream(new File(outFileName));
            } catch (IOException e) {
                System.err.println("Error opening file: "+outFileName);
                System.exit(-1);
            }
        }

        if (!summaryOnly) {
            outputHeader(traits);

            if (slices == null)
                outputSlice(0,Double.NaN);
            else {
                for(int i=0; i<slices.length; i++)
                    outputSlice(i,slices[i]);

            }
        } else { // Output summaries

            if (outputFormat == OutputFormat.XML) {
                rootElement = new Element("xml");

            } else if (outputFormat == OutputFormat.KML) {
                Element folderNameElement = new Element("name");
                folderNameElement.addContent("surface HPD regions");

                folderElement = new Element("Folder");
                //rootElement.setAttribute("xmlns","http://earth.google.com/kml/2.2");
                folderElement.addContent(folderNameElement);

                Element documentNameElement = new Element("name");
                String documentName = outFileName;
                if (documentName == null)
                    documentName = "default";
                if (documentName.endsWith(".kml"))
                    documentName = documentName.replace(".kml","");
                documentNameElement.addContent(documentName);

                documentElement = new Element("Document");
                documentElement.addContent(documentNameElement);
                documentElement.addContent(folderElement);

                rootElement = new Element("kml");
                rootElement.addContent(documentElement);
            }

            if (slices == null)
                summarizeSlice(0,Double.NaN, outputFormat, hpdValue);
            else {
                for(int i=0; i<slices.length; i++)
                    summarizeSlice(i,slices[i], outputFormat, hpdValue);
            }

            XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
            try {
                xmlOutputter.output(rootElement,resultsStream);
            } catch (IOException e) {
                System.err.println("IO Exception encountered: "+e.getMessage());
                System.exit(-1);
            }
        }
// writes out dispersal rate summaries for the whole tree, there is a beast xml statistic that can do this now
//        if (containsLocation) {
//            try{
//                PrintWriter dispersalRateFile = new PrintWriter(new FileWriter("dispersalRates.log"), true);
//                dispersalRateFile.print("state"+"\t"+"dispersalRate(native units)"+"\t"+"dispersalRate(great circle distance, km)\r");
//                for(int x = 0; x < dispersalrates.size(); x++ ) {
//                    dispersalRateFile.print(x+"\t"+dispersalrates.get(x)+"\r");
//                }
//
//            } catch (IOException e) {
//                System.err.println("IO Exception encountered: "+e.getMessage());
//                System.exit(-1);
//            }
//        }


//  <--          attempt to get slice dispersalRates
            if (sdr) {
                double[][] sliceTreeDistances = new double[sliceCount][sliceTreeDistanceArrays.size()];
                double[][] sliceTreeTimes = new double[sliceCount][sliceTreeDistanceArrays.size()];
                double[][] sliceTreeRates = new double[sliceTreeDistanceArrays.size()][sliceCount];
                for (int q = 0; q < sliceTreeDistanceArrays.size(); q++){
                    double[] distanceArray = (double[])sliceTreeDistanceArrays.get(q);
                    double[] timeArray = (double[])sliceTreeTimeArrays.get(q);
                    for (int r = 0; r < distanceArray.length; r ++) {
                        sliceTreeDistances[r][q] = distanceArray[r];
                        sliceTreeTimes[r][q] = timeArray[r];
                    }
                }

                print2DArray(sliceTreeDistances,"sliceTreeDistances.txt");
                print2DArray(sliceTreeTimes,"sliceTreeTimes.txt");

                if (sliceCount > 1) {
                    for (int s = 0; s < sliceTreeDistanceArrays.size(); s++) {
                        double[] distanceArray = (double[])sliceTreeDistanceArrays.get(s);
                        double[] timeArray = (double[])sliceTreeTimeArrays.get(s);
                        for (int t = 0; t <  (sliceCount-1); t++) {
                            sliceTreeRates[s][t] = (distanceArray[t] - distanceArray[t+1])/(timeArray[t] - timeArray[t+1]);
                            //sliceTreeRates[s][t] = (sliceTreeDistances[t][s] - sliceTreeDistances[t+1][s])/(sliceTreeTimes[t][s] - sliceTreeTimes[t+1][s]);
                        }
                        sliceTreeRates[s][sliceCount-1] = (sliceTreeDistances[sliceCount-1][s])/(sliceTreeTimes[sliceCount-1][s]);
                    }
                } else {
                    for (int s = 0; s < sliceTreeDistanceArrays.size(); s++) {
                        sliceTreeRates[s][0] = sliceTreeDistances[0][s]/sliceTreeTimes[0][s];
                    }
                }

                print2DArray(sliceTreeRates,"sliceTreeRates.txt");

                try{
                    PrintWriter sliceDispersalRateFile = new PrintWriter(new FileWriter(sdrFile), true);
                    sliceDispersalRateFile.print("sliceTime"+"\t");
                    if (mostRecentSamplingDate > 0) {
                        sliceDispersalRateFile.print("realTime"+"\t");
                    }
                    sliceDispersalRateFile.print("mean dispersalRate"+"\t"+"\t"+"hpd low"+"\t"+"hpd up"+"\r");
                    double[] meanDispersalRates = meanColNoNaN(sliceTreeRates);
                    double[][] hpdDispersalRates = getArrayHPDintervals(sliceTreeRates);
                    for (int u = 0; u < sliceCount; u++) {
                        sliceDispersalRateFile.print(slices[u]+"\t");
                        if (mostRecentSamplingDate > 0) {
                            sliceDispersalRateFile.print((mostRecentSamplingDate-slices[u])+"\t");
                        }
                        sliceDispersalRateFile.print(meanDispersalRates[u] + "\t"+ hpdDispersalRates[u][0] + "\t" + hpdDispersalRates[u][1]+"\r");
                    }
                    sliceDispersalRateFile.close();
                } catch (IOException e) {
                    System.err.println("IO Exception encountered: "+e.getMessage());
                    System.exit(-1);
                }
//               attempt to get slice dispersalRates  -->
            }
    }

    enum Normalization {
        LENGTH,
        HEIGHT,
        NONE
    }

    enum OutputFormat {
        TAB,
        KML,
        XML
    }

    public static <T extends Enum<T>> String[] enumNamesToStringArray(T[] values) {
        int i = 0;
        String[] result = new String[values.length];
        for (T value: values) {
            result[i++] = value.name();
        }
        return result;
    }



    private void addDimInfo(Element element, int j, int dim) {
        if (dim > 1)
            element.setAttribute("dim",Integer.toString(j+1));
    }

    private void summarizeSliceTrait(Element sliceElement, int slice, List<Trait> thisTrait, int traitIndex, double sliceValue,
                                     OutputFormat outputFormat,
                                     double hpdValue) {

                boolean isNumber = thisTrait.get(0).isNumber();
                boolean isMultivariate = thisTrait.get(0).isMultivariate();
                int dim = thisTrait.get(0).getValue().length;
                boolean isBivariate = isMultivariate && dim == 2;
                if (isNumber) {

                    Element traitElement = null;
                    if (outputFormat == OutputFormat.XML) {
                        traitElement = new Element(TRAIT);
                        traitElement.setAttribute(NAME,traits[traitIndex]);
                    }

                    if (outputFormat == OutputFormat.KML) {
                        Element styleElement = new Element(STYLE);
                        constructPolygonStyleElement(styleElement, sliceValue);
                        documentElement.addContent(styleElement);
                    }

                    int count = thisTrait.size();
                    double[][] x = new double[dim][count];
                    for(int i=0; i<count; i++) {
                        Trait trait = thisTrait.get(i);
                        double[] value = trait.getValue();
                        for(int j=0; j<dim; j++)
                            x[j][i] = value[j];
                    }

                    if (outputFormat == OutputFormat.XML) {
                        // Compute marginal means and standard deviations
                        for (int j = 0; j < dim; j++) {
                            TraceDistribution trace = new TraceDistribution(x[j]);
                            Element statsElement = new Element("stats");
                            addDimInfo(statsElement, j, dim);
                            StringBuffer sb = new StringBuffer();
                            sb.append(KMLCoordinates.NEWLINE);
                            sb.append(String.format(KMLCoordinates.FORMAT,
                                    trace.getMean())).append(KMLCoordinates.SEPERATOR);
                            sb.append(String.format(KMLCoordinates.FORMAT,
                                    trace.getStdError())).append(KMLCoordinates.SEPERATOR);
                            sb.append(String.format(KMLCoordinates.FORMAT,
                                    trace.getLowerHPD())).append(KMLCoordinates.SEPERATOR);
                            sb.append(String.format(KMLCoordinates.FORMAT,
                                    trace.getUpperHPD())).append(KMLCoordinates.NEWLINE);
                            statsElement.addContent(sb.toString());
                            traitElement.addContent(statsElement);
                        }
                    }

                    if (isBivariate) {

                        //for testing how much points are within the polygons
                        double numberOfPointsInPolygons = 0;
                        double totalArea = 0;

                        ContourMaker contourMaker;
                        if (contourMode == ContourMode.JAVA)
                            contourMaker = new KernelDensityEstimator2D(x[0],x[1]);
                        else if (contourMode == ContourMode.R)
                            contourMaker = new ContourWithR(x[0],x[1]);
                        else if (contourMode == ContourMode.SNYDER)
                            contourMaker = new ContourWithSynder(x[0],x[1]);
                        else
                            throw new RuntimeException("Unimplemented ContourModel!");

                        ContourPath[] paths = contourMaker.getContourPaths(hpdValue);
                        for(ContourPath path : paths) {

                            KMLCoordinates coords = new KMLCoordinates(path.getAllX(),path.getAllY());
                            if(outputFormat == OutputFormat.XML) {
                                Element regionElement = new Element(REGIONS_ELEMENT);
                                regionElement.setAttribute(DENSITY_VALUE,Double.toString(hpdValue));
                                regionElement.addContent(coords.toXML());
                                traitElement.addContent(regionElement);
                            }

                            // only if the trait is location we will write KML
                            if(outputFormat == OutputFormat.KML){
                                //because KML polygons require long,lat,alt we need to switch lat and long first
                                coords.switchXY();
                                Element placemarkElement = generatePlacemarkElementWithPolygon(sliceValue, coords, slice);
                                //testing how many points are within the polygon
                                Element testElement = new Element("test");
                                testElement.addContent(coords.toXML());
                                Polygon2D testPolygon = new Polygon2D(testElement);
                                totalArea += testPolygon.calculateArea();
                                numberOfPointsInPolygons += getNumberOfPointsInPolygon(x,testPolygon);

                                folderElement.addContent(placemarkElement);
                            }
                       }
                       //testing how many points are within the polygon
                       System.out.println(sliceValue+"\t"+(mostRecentSamplingDate-sliceValue)+"\t"+paths.length+"\t"+numberOfPointsInPolygons/count+"\t"+totalArea);

                    }
                    if (outputFormat == OutputFormat.XML)
                        sliceElement.addContent(traitElement);

                } // else skip

    }

    public static int getNumberOfPointsInPolygon(double[][] pointsArray, Polygon2D testPolygon) {
        int numberOfPointsInPolygon = 0;
        for (int x = 0; x < pointsArray[0].length; x++){
            if (testPolygon.containsPoint2D(new Point2D.Double(pointsArray[1][x],pointsArray[0][x]))){
                numberOfPointsInPolygon++;
            }
        }
        return numberOfPointsInPolygon;
    }

    private void summarizeSlice(int slice, double sliceValue, OutputFormat outputFormat, double hpdValue) {

        if (outputFormat == OutputFormat.TAB)
            throw new RuntimeException("Only XML/KML output is implemented");

            Element sliceElement = null;

            if (outputFormat == OutputFormat.XML) {
                sliceElement = new Element(SLICE_ELEMENT);
                sliceElement.setAttribute(SLICE_VALUE, Double.toString(sliceValue));
            }

            List<List<Trait>> thisSlice = values.get(slice);
            int traitCount = thisSlice.size();

            for(int traitIndex=0; traitIndex<traitCount; traitIndex++) {

                summarizeSliceTrait(sliceElement, slice, thisSlice.get(traitIndex), traitIndex, sliceValue,
                                     outputFormat,
                                     hpdValue);
            }

            if(outputFormat == OutputFormat.XML) {
                rootElement.addContent(sliceElement);
            }
    }

    private void constructPolygonStyleElement(Element styleElement, double sliceValue){
        double date = mostRecentSamplingDate - sliceValue;
        styleElement.setAttribute(ID,REGIONS_ELEMENT+date+"_style");
            Element lineStyle = new Element("LineStyle");
                Element width = new Element("width");
                width.addContent(WIDTH);
            lineStyle.addContent(width);
            Element polyStyle = new Element("PolyStyle");
                Element color = new Element("color");
                double[] minMax = new double[2];
                minMax[0] = slices[0];
                minMax[1] = slices[(slices.length-1)];
                String colorString = getKMLColor(sliceValue,minMax,endHPDColor,startHPDColor);
                color.addContent(opacity+colorString);
                Element outline = new Element("outline");
                outline.addContent("0");
            polyStyle.addContent(color);
            polyStyle.addContent(outline);
        styleElement.addContent(lineStyle);
        styleElement.addContent(polyStyle);
    }

    private Element generatePlacemarkElementWithPolygon(double sliceValue, KMLCoordinates coords, int sliceInteger){
        double date = mostRecentSamplingDate - sliceValue;
        Element placemarkElement  = new Element("Placemark");

        Element visibility = new Element("visibility");
        visibility.addContent("1");
        placemarkElement.addContent(visibility);

        if (sliceCount > 1) {
            Element timeSpan = new Element("TimeSpan");
                Element begin = new Element("begin");
                begin.addContent(Double.toString(date));
                timeSpan.addContent(begin);
            if (sliceInteger > 1) {
                 Element end = new Element("end");
                end.addContent(Double.toString(mostRecentSamplingDate-slices[(sliceInteger-1)]));
                timeSpan.addContent(end);
            }
            placemarkElement.addContent(timeSpan);
        }

        Element style = new Element("styleUrl");
        style.addContent("#"+REGIONS_ELEMENT+date+"_style");
        placemarkElement.addContent(style);

        Element polygonElement  = new Element("Polygon");
        Element altitudeMode = new Element("altitudeMode");
        altitudeMode.addContent("clampToGround");
        polygonElement.addContent(altitudeMode);
        Element tessellate = new Element("tessellate");
        tessellate.addContent("1");
        polygonElement.addContent(tessellate);
        Element outerBoundaryIs = new Element("outerBoundaryIs");
        Element LinearRing = new Element("LinearRing");
        LinearRing.addContent(coords.toXML());
        outerBoundaryIs.addContent(LinearRing);
        polygonElement.addContent(outerBoundaryIs);
        placemarkElement.addContent(polygonElement);

        return placemarkElement;
    }

    private void outputHeader(String[] traits) {
        StringBuffer sb = new StringBuffer("slice");
        for(int i=0; i<traits.length; i++) {
            // Load first value to check dimensionality
            Trait trait = values.get(0).get(i).get(0);
            if (trait.isMultivariate()) {
                int dim = trait.getDim();
                for(int j=1; j<=dim; j++)
                    sb.append(sep).append(traits[i]).append(j);
            } else
                sb.append(sep).append(traits[i]);
        }
        sb.append("\n");
        resultsStream.print(sb);
    }

//    private List<Tree> importTrees(String treeFileName, int burnin) throws IOException, Importer.ImportException {
//
//        int totalTrees = 10000;
//
//        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
//        progressStream.println("0              25             50             75            100");
//        progressStream.println("|--------------|--------------|--------------|--------------|");
//
//        int stepSize = totalTrees / 60;
//        if (stepSize < 1) stepSize = 1;
//
//        List<Tree> treeList = new ArrayList<Tree>();
//        BufferedReader reader1 = new BufferedReader(new FileReader(treeFileName));
//
//        String line1 = reader1.readLine();
//        TreeImporter importer1;
//        if (line1.toUpperCase().startsWith("#NEXUS")) {
//            importer1 = new NexusImporter(new FileReader(treeFileName));
//        } else {
//            importer1 = new NewickImporter(new FileReader(treeFileName));
//        }
//        totalTrees = 0;
//        while (importer1.hasTree()) {
//            Tree treeTime = importer1.importNextTree();
//
//            if (totalTrees > burnin)
//                treeList.add(treeTime);
//
//            if (totalTrees > 0 && totalTrees % stepSize == 0) {
//                progressStream.print("*");
//                progressStream.flush();
//            }
//            totalTrees++;
//        }
//        return treeList;
//    }

    private void readAndAnalyzeTrees(String treeFileName, int burnin,
                                     String[] traits, double[] slices,
                                     boolean impute, boolean trueNoise, Normalization normalize)
            throws IOException, Importer.ImportException {

        int totalTrees = 10000;
        int totalStars = 0;

        progressStream.println("Reading and analyzing trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        BufferedReader reader1 = new BufferedReader(new FileReader(treeFileName));

        String line1 = reader1.readLine();
        TreeImporter importer1;
        if (line1.toUpperCase().startsWith("#NEXUS")) {
            importer1 = new NexusImporter(new FileReader(treeFileName));
        } else {
            importer1 = new NewickImporter(new FileReader(treeFileName));
        }
        totalTrees = 0;
        while (importer1.hasTree()) {
            Tree treeTime = importer1.importNextTree();
            treesRead++;
            if (totalTrees > burnin)
               analyzeTree(treeTime, traits, slices, impute, trueNoise, normalize);

            if (totalTrees > 0 && totalTrees % stepSize == 0) {
                progressStream.print("*");
                totalStars++;
                if(totalStars % 61 == 0)
                    progressStream.print("\n");
                progressStream.flush();
            }
            totalTrees++;
        }
        progressStream.print("\n");        
    }

    class Trait {

        Trait(Object obj) {
            this.obj = obj;
            if (obj instanceof Object[]) {
                isMultivariate = true;
                array = (Object[])obj;
            }
        }

        public boolean isMultivariate() { return isMultivariate; }

        public boolean isNumber() {
            if (!isMultivariate)
                return (obj instanceof Double);
            return (array[0] instanceof Double);
        }

        public int getDim() {
            if (isMultivariate) {
                return array.length;
            }
            return 1;
        }

        public double[] getValue() {
            int dim = getDim();
            double[] result = new double[dim];
            for(int i=0; i<dim; i++)
                result[i] = (Double)array[i];
            return result;
        }

        private Object obj;
        private Object[] array;
        private boolean isMultivariate = false;

        public String toString() {
            if (!isMultivariate)
                return obj.toString();
            StringBuffer sb = new StringBuffer(array[0].toString());
            for(int i=1; i<array.length; i++)
                sb.append(sep).append(array[i]);
            return sb.toString();
        }
    }

    private List<List<List<Trait>>> values;

    private void outputSlice(int slice, double sliceValue) {

        List<List<Trait>> thisSlice = values.get(slice);
        int traitCount = thisSlice.size();
        int valueCount = thisSlice.get(0).size();

        StringBuffer sb = new StringBuffer();

        for(int v=0; v<valueCount; v++) {
            if (Double.isNaN(sliceValue))
                sb.append("All");
            else
                sb.append(sliceValue);
            for(int t=0; t<traitCount; t++) {
                sb.append(sep);
                sb.append(thisSlice.get(t).get(v));
            }
            sb.append("\n");
        }

       resultsStream.print(sb);            
    }

    private void analyzeTree(Tree treeTime, String[] traits, double[] slices, boolean impute,
                             boolean trueNoise, Normalization normalize) {

        double[][] precision = null;

        if (impute) {

            Object o = treeTime.getAttribute(PRECISION_STRING);
            double treeNormalization = 1; // None
            if (normalize == Normalization.LENGTH) {
                treeNormalization = Tree.Utils.getTreeLength(treeTime, treeTime.getRoot());
            } else if (normalize == Normalization.HEIGHT) {
                treeNormalization = treeTime.getNodeHeight(treeTime.getRoot());
            }

            if (o != null) {
                Object[] array = (Object[]) o;
                int dim = (int) Math.sqrt(1 + 8 * array.length) / 2;
                precision = new double[dim][dim];
                int c = 0;
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        precision[j][i] = precision[i][j] = ((Double) array[c++])*treeNormalization;
                    }
                }
            }
        }

//  employed to get dispersal rates across the whole tree
//        double treeNativeDistance = 0;
//        double treeKilometerGreatCircleDistance = 0;

//  <--  attempt to get slice dispersalRates
        double[] treeSliceTime = new double[sliceCount];
        double[] treeSliceDistance = new double[sliceCount];
//       attempt to get slice dispersalRates  -->

        for (int x = 0; x < treeTime.getNodeCount(); x++) {

            NodeRef node = treeTime.getNode(x);


            if (!(treeTime.isRoot(node))) {

                double nodeHeight = treeTime.getNodeHeight(node);
                double parentHeight = treeTime.getNodeHeight(treeTime.getParent(node));

//  employed to get dispersal rates across the whole tree
//                if (containsLocation){
//
//                    Trait nodeLocationTrait = new Trait (treeTime.getNodeAttribute(node, LOCATIONTRAIT));
//                    Trait parentNodeLocationTrait = new Trait (treeTime.getNodeAttribute(treeTime.getParent(node), LOCATIONTRAIT));
//                    treeNativeDistance += getNativeDistance(nodeLocationTrait.getValue(),parentNodeLocationTrait.getValue());
//                    treeKilometerGreatCircleDistance += getKilometerGreatCircleDistance(nodeLocationTrait.getValue(),parentNodeLocationTrait.getValue());
//
//                }

                for (int i = 0; i < sliceCount; i++) {

//  <--          attempt to get slice dispersalRates
                    if (sdr) {
                        if (!doSlices ||
                                (slices[i] <= nodeHeight)
                                ) {
                                treeSliceTime[i] += (parentHeight-nodeHeight);
                                //TreeModel model = new TreeModel(treeTime, true);
                                //treeSliceDistance[i] += getKilometerGreatCircleDistance(model.getMultivariateNodeTrait(node, LOCATIONTRAIT),model.getMultivariateNodeTrait(model.getParent(node), LOCATIONTRAIT));
                                Trait nodeLocationTrait = new Trait (treeTime.getNodeAttribute(node, LOCATIONTRAIT));
                                Trait parentNodeLocationTrait = new Trait (treeTime.getNodeAttribute(treeTime.getParent(node), LOCATIONTRAIT));
                                treeSliceDistance[i] += getKilometerGreatCircleDistance(nodeLocationTrait.getValue(),parentNodeLocationTrait.getValue());
                        }
                    }
//       attempt to get slice dispersalRates  -->

                    if (!doSlices ||
                            (slices[i] >= nodeHeight && slices[i] < parentHeight)
                            ) {

                        List<List<Trait>> thisSlice = values.get(i);
                        for (int j = 0; j < traitCount; j++) {

                            List<Trait> thisTraitSlice = thisSlice.get(j);
                            Object tmpTrait = treeTime.getNodeAttribute(node, traits[j]);
                            if (tmpTrait == null) {
                                System.err.println("Trait '" + traits[j] + "' not found on branch.");
                                System.exit(-1);
                            }
                            Trait trait = new Trait(tmpTrait);
                            if (impute) {
                                Double rateAttribute = (Double) treeTime.getNodeAttribute(node, RATE_STRING);
                                double rate = 1.0;
                                if (rateAttribute != null) {
                                    rate = rateAttribute;
                                    if (outputRateWarning) {
                                        progressStream.println("Warning: using a rate attribute during imputation!");
                                        outputRateWarning = false;
                                    }
                                }
                                if (trueNoise && precision == null) {
                                    progressStream.println("Error: not precision available for imputation with correct noise!");
                                    System.exit(-1);
                                }
                                trait = imputeValue(trait, new Trait(treeTime.getNodeAttribute(treeTime.getParent(node), traits[j])),
                                        slices[i], nodeHeight, parentHeight, precision, rate, trueNoise);
                            }
                            thisTraitSlice.add(trait);

//  <--          attempt to get slice dispersalRates
                 //if trait is location
                            if (sdr){
                                treeSliceTime[i] += (parentHeight-slices[i]);
                                Trait parentTrait = new Trait(treeTime.getNodeAttribute(treeTime.getParent(node), traits[j]));
                                treeSliceDistance[i] += getKilometerGreatCircleDistance(trait.getValue(), parentTrait.getValue());
//               attempt to get slice dispersalRates  -->
                            }
                        }
                    }
                }
            }
        }

//  <--          attempt to get slice dispersalRates
        if (sdr){
            sliceTreeDistanceArrays.add(treeSliceDistance);
            sliceTreeTimeArrays.add(treeSliceTime);
        }
//               attempt to get slice dispersalRates  -->
        

//  employed to get dispersal rates across the whole tree
//        if (containsLocation) {
//           double treelength = Tree.Utils.getTreeLength(treeTime, treeTime.getRoot());
//            double dispersalNativeRate = treeNativeDistance/treelength;
//            double dispersalKilometerRate = treeKilometerGreatCircleDistance/treelength;
            //System.out.println(dispersalNativeRate+"\t"+dispersalKilometerRate);
//            dispersalrates.add(dispersalNativeRate+"\t"+dispersalKilometerRate);
//        }

        treesAnalyzed++;

    }
//  employed to get dispersal rates across the whole tree
//    private static double getNativeDistance(double[] location1, double[] location2) {
//        return Math.sqrt(Math.pow((location2[0]-location1[0]),2.0)+Math.pow((location2[1]-location1[1]),2.0));
//    }
//
   private static double getKilometerGreatCircleDistance(double[] location1, double[] location2) {
        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(location1[0], location1[1]);
        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(location2[0], location2[1]);
        return (coord1.distance(coord2));
   }

    private int traitCount;
    private int sliceCount;
    private String[] traits;
    private double[] slices;
    private boolean doSlices;
    private int treesRead = 0;
    private int treesAnalyzed = 0;
    private double mostRecentSamplingDate;
    private ContourMode contourMode;
//  employed to get dispersal rates across the whole tree
//    private static boolean containsLocation = false;
//    private static ArrayList dispersalrates = new ArrayList();

//    private void run(List<Tree> trees, String[] traits, double[] slices, boolean impute, boolean trueNoise) {
//
//        for (Tree treeTime : trees) {
//            analyzeTree(treeTime, traits, slices, impute, trueNoise);
//        }
//    }

//  <--  attempt to get slice dispersalRates
    private ArrayList sliceTreeDistanceArrays = new ArrayList();
    private ArrayList sliceTreeTimeArrays = new ArrayList();
//    attempt to get slice dispersalRates  -->
    private boolean sdr;

    private boolean outputRateWarning = true;


    private Trait imputeValue(Trait nodeTrait, Trait parentTrait, double time, double nodeHeight, double parentHeight, double[][] precision, double rate, boolean trueNoise) {
        if (!nodeTrait.isNumber()) {
            System.err.println("Can only impute numbers!");
            System.exit(-1);
        }

        int dim = nodeTrait.getDim();
        double[] nodeValue = nodeTrait.getValue();
        double[] parentValue = parentTrait.getValue();

        final double timeTotal = parentHeight - nodeHeight;
        final double timeChild = (time - nodeHeight);
        final double timeParent = (parentHeight - time);
        final double weightTotal = 1.0 / timeChild + 1.0 / timeParent;

        if (timeChild == 0)
            return nodeTrait;

        if (timeParent == 0)
            return parentTrait;

        // Find mean value, weighted average
        double[] mean = new double[dim];
        double[][] scaledPrecision = new double[dim][dim];

        for(int i=0; i<dim; i++) {
            mean[i] = (nodeValue[i] / timeChild + parentValue[i] / timeParent) / weightTotal;
            if (trueNoise) {
                for(int j=i; j<dim; j++)
                    scaledPrecision[j][i] = scaledPrecision[i][j] = precision[i][j] / timeTotal / rate;
            }
        }

        if (trueNoise)
            mean = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precision);
        Object[] result = new Object[dim];
        for(int i=0; i<dim; i++)
            result[i] = mean[i];
        return new Trait(result);
    }

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;
    private PrintStream resultsStream;

    private final static Version version = new BeastVersion();

    private static final String commandName = "timeslicer";


    public static void printUsage(Arguments arguments) {

        arguments.printUsage(commandName, "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: " + commandName + " test.trees out.kml");
        progressStream.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            progressStream.print(" ");
        }
        progressStream.println(line);
    }

    public static void printTitle() {
        progressStream.println();
        centreLine("TimeSlicer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Marc A. Suchard, Philippe Lemey,", 60);
        centreLine("Alexei J. Drummond and Andrew Rambaut", 60);
        progressStream.println();
        centreLine("Department of Biomathematics", 60);
        centreLine("University of California, Los Angeles", 60);
        centreLine("msuchard@ucla.edu", 60);
        progressStream.println();
        centreLine("Rega Institute for Medical Research", 60);
        centreLine("Katholieke Universiteit Leuven", 60);
        centreLine("philippe.lemey@gmail.com", 60);
        progressStream.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        progressStream.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        progressStream.println();
        progressStream.println();
    }

    private static double[] parseVariableLengthDoubleArray(String inString) throws Arguments.ArgumentException {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString,",");
        while(st.hasMoreTokens()) {
            try {
                returnList.add(Double.parseDouble(st.nextToken()));
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }

        }

        if (returnList.size()>0) {
            double[] doubleArray = new double[returnList.size()];
            for(int i=0; i<doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);
            return doubleArray;
        }
        return null;
    }

    private static String[] parseVariableLengthStringArray(String inString) {

       List<String> returnList = new ArrayList<String>();
       StringTokenizer st = new StringTokenizer(inString,",");
       while(st.hasMoreTokens()) {
           returnList.add(st.nextToken());
       }

       if (returnList.size()>0) {
           String[] stringArray = new String[returnList.size()];
           stringArray = returnList.toArray(stringArray);
           return stringArray;
       }
       return null;
   }

    private static double[] parseFileWithArray(String file) {
        List<Double> returnList = new ArrayList<Double>();
        try{
            BufferedReader readerTimes = new BufferedReader(new FileReader(file));
            String line = readerTimes.readLine();
            while (line != null && !line.equals("")) {
                returnList.add(Double.valueOf(line));
                line = readerTimes.readLine();
            }
        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        if (returnList.size()>0) {
            double[] doubleArray = new double[returnList.size()];
            for(int i=0; i<doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);
            return doubleArray;
        }
        return null;
   }

    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String[] traitNames = null;
        double[] sliceTimes = null;
        OutputFormat outputFormat = OutputFormat.XML;
        boolean impute = false;
        boolean trueNoise = true;
        boolean summaryOnly = false;
        ContourMode contourMode = ContourMode.SNYDER;
        Normalization normalize = Normalization.LENGTH;
        int burnin = -1;
        double mrsd = 0;
        double hpdValue = 0.80;
        String outputFileSDR = null;
        boolean getSDR = false;


//        if (args.length == 0) {
//          // TODO Make flash GUI
//        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURNIN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption(TRAIT, "trait_name", "specifies an attribute-list to use to create a density map [default = location.rate]"),
                        new Arguments.StringOption(SLICE,"time","specifies an slice time-list [default=none]"),
                        new Arguments.StringOption(FILE_SLICE,"fileTimes","specifies a file with a slice time-list, is overwritten by command-line specification of slice times [default=none]"),
                        new Arguments.RealOption(MRSD,"specifies the most recent sampling data in fractional years to rescale time [default=0]"),
                        new Arguments.Option(HELP, "option to print this message"),
                        new Arguments.StringOption(NOISE, falseTrue, false,
                                "add true noise [default = true])"),
                        new Arguments.StringOption(IMPUTE, falseTrue, false,
                                "impute trait at time-slice [default = false]"),
                        new Arguments.StringOption(SUMMARY, falseTrue, false,
                                "compute summary statistics [default = false]"),
                        new Arguments.StringOption(FORMAT, enumNamesToStringArray(OutputFormat.values()),false,
                                "summary output format [default = XML]"),
                        new Arguments.IntegerOption(HPD,"mass (1 - 99%) to include in HPD regions [default = 80]"),
                        new Arguments.StringOption(MODE,enumNamesToStringArray(ContourMode.values()), false,
                                "contouring model [default = synder]"),
                        new Arguments.StringOption(NORMALIZATION, enumNamesToStringArray(Normalization.values()), false,
                                "tree normalization [default = length"),
                        new Arguments.StringOption(SDR,"sliceDispersalRate","specifies output file name for dispersal rates for each slice (from previous sliceTime[or root of the trees] up to current sliceTime"),

                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption(HELP)) {
            printUsage(arguments);
            System.exit(0);
        }

        try { // Make sense of arguments

            String sliceFileString = arguments.getStringOption(FILE_SLICE);
            if (sliceFileString != null) {
                sliceTimes = parseFileWithArray(sliceFileString);
            }


            if (arguments.hasOption(BURNIN)) {
                burnin = arguments.getIntegerOption(BURNIN);
            }


            if (arguments.hasOption(MRSD)) {
                mrsd = arguments.getRealOption(MRSD);
            }


            if (arguments.hasOption(HPD)) {
                int intValue = arguments.getIntegerOption(HPD);
                if (intValue < 1 || intValue > 99) {
                    progressStream.println("HPD Region mass falls outside of 1 - 99% range.");
                    System.exit(-1);
                }
                hpdValue = intValue / 100.0;
            }

            String traitString = arguments.getStringOption(TRAIT);
            if (traitString != null) {
                traitNames = parseVariableLengthStringArray(traitString);
//employed to get dispersal rates across the whole tree
//                for (int y = 0; y < traitNames.length; y++) {
//                    if (traitNames[y].equals(LOCATIONTRAIT)){
//                       containsLocation =  true;
//                    }
//                }
            }

            if (traitNames == null) {
                traitNames = new String[1];
                traitNames[0] = "location.rate";
            }

            String sliceString = arguments.getStringOption(SLICE);
            if (sliceString != null) {
                sliceTimes = parseVariableLengthDoubleArray(sliceString);
            }

            String imputeString = arguments.getStringOption(IMPUTE);
            if (imputeString != null && imputeString.compareToIgnoreCase("true") == 0)
                impute = true;

            String noiseString = arguments.getStringOption(NOISE);
            if (noiseString != null && noiseString.compareToIgnoreCase("false") == 0)
                trueNoise = false;

            String summaryString = arguments.getStringOption(SUMMARY);
            if (summaryString != null && summaryString.compareToIgnoreCase("true") == 0)
                summaryOnly = true;

            String modeString = arguments.getStringOption(MODE);
            if (modeString != null) {
                contourMode = ContourMode.valueOf(modeString.toUpperCase());
                if (contourMode == ContourMode.R && !ContourWithR.processWithR)
                    contourMode = ContourMode.SNYDER;
            }

            String normalizeString = arguments.getStringOption(NORMALIZATION);
            if (normalizeString != null) {
                normalize = Normalization.valueOf(normalizeString.toUpperCase());
            }

            String summaryFormat = arguments.getStringOption(FORMAT);
            if (summaryFormat != null) {
                outputFormat = OutputFormat.valueOf(summaryFormat.toUpperCase());
            }

            String sdrString = arguments.getStringOption("sdr");
            if (sdrString != null) {
                outputFileSDR = sdrString;
                getSDR = true;
            }

        } catch (Arguments.ArgumentException e) {
            progressStream.println(e);
            printUsage(arguments);
            System.exit(-1);
        }


        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 0:
                printUsage(arguments);
                System.exit(1);
            case 2:
                outputFileName = args2[1];
                // fall to
            case 1:
                inputFileName = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        TimeSlicer timeSlicer = new TimeSlicer(inputFileName, burnin, traitNames, sliceTimes, impute,
                trueNoise, mrsd, contourMode, normalize, getSDR);
        timeSlicer.output(outputFileName, summaryOnly, outputFormat, hpdValue, outputFileSDR);

        System.exit(0);
    }




    public static String getKMLColor(double value, double[] minMaxMedian, String startColor, String endColor) {

        startColor = startColor.toLowerCase();
        String startBlue = startColor.substring(0,2);
        String startGreen = startColor.substring(2,4);
        String startRed = startColor.substring(4,6);

        endColor =  endColor.toLowerCase();
        String endBlue = endColor.substring(0,2);
        String endGreen = endColor.substring(2,4);
        String endRed = endColor.substring(4,6);

        double proportion = (value - minMaxMedian[0])/(minMaxMedian[1] - minMaxMedian[0]);

        // generate an array with hexadecimal code for each RGB entry number
        String[] colorTable = new String[256];

        int colorTableCounter = 0;

        for (int a = 0; a < 10; a++) {

            for (int b = 0; b < 10; b++) {

                colorTable[colorTableCounter] = a + "" + b;
                colorTableCounter ++;
            }

            for(int c = (int)('a'); c<6+(int)('a'); c++) {
                colorTable[colorTableCounter] = a + "" + (char)c;
                colorTableCounter ++;
            }

        }
        for(int d = (int)('a'); d<6+(int)('a'); d++) {

            for (int e = 0; e < 10; e++) {

                colorTable[colorTableCounter] = (char) d + "" + e;
                colorTableCounter ++;
            }

            for(int f = (int)('a'); f<6+(int)('a'); f++) {
                colorTable[colorTableCounter] = (char) d + "" + (char) f;
                colorTableCounter ++;
            }

        }


        int startBlueInt = 0;
        int startGreenInt = 0;
        int startRedInt = 0;

        int endBlueInt = 0;
        int endGreenInt = 0;
        int endRedInt = 0;

        for (int i = 0; i < colorTable.length; i ++) {

            if (colorTable[i].equals(startBlue)) {startBlueInt = i; }
            if (colorTable[i].equals(startGreen)) {startGreenInt = i; }
            if (colorTable[i].equals(startRed)) {startRedInt = i; }
            if (colorTable[i].equals(endBlue)) {endBlueInt = i; }
            if (colorTable[i].equals(endGreen)) {endGreenInt = i; }
            if (colorTable[i].equals(endRed)) {endRedInt = i; }

        }

        int blueInt = startBlueInt + (int) Math.round((endBlueInt-startBlueInt)*proportion);
        int greenInt = startGreenInt + (int) Math.round((endGreenInt-startGreenInt)*proportion);
        int redInt = startRedInt + (int) Math.round((endRedInt-startRedInt)*proportion);

        String blue = null;
        String green =  null;
        String red = null;

        for (int j = 0; j < colorTable.length; j ++) {

            if (j == blueInt) {blue = colorTable[j]; }
            if (j == greenInt) {green = colorTable[j]; }
            if (j == redInt) {red = colorTable[j]; }

        }

        return blue+green+red;
    }

    private static double[] getHPDInterval(double proportion, double[] array, int[] indices) {

        double returnArray[] = new double[2];
        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int)Math.round(proportion * (double)array.length);
        for (int i = 0; i <= (array.length - diff); i++) {
            double minValue = array[indices[i]];
            double maxValue = array[indices[i+diff-1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        returnArray[0] = array[indices[hpdIndex]];
        returnArray[1] = array[indices[hpdIndex+diff-1]];
        return returnArray;
    }

    private static void print2DArray(double[][] array, String name) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[0].length; j++) {
                    outFile.print(array[i][j]+"\t");
                }
            outFile.println("");
            }
            outFile.close();

        } catch(IOException io) {
           System.err.print("Error writing to file: " + name);
        }
    }

    private static double[] meanColNoNaN(double[][] x)    {
        double[] returnArray = new double[x[0].length];

        for (int i = 0; i < x[0].length; i++) {

         	double m = 0.0;
        	int lenNoZero = 0;


        	for (int j = 0; j < x.length; j++) {

                if (!(((Double)x[j][i]).isNaN())) {
                    m += x[j][i];
                    lenNoZero += 1;
                }
        	}
        	returnArray[i] = m / (double) lenNoZero;

        }
        return returnArray;
    }

    private static double[][] getArrayHPDintervals(double[][] array) {

        double[][] returnArray = new double[array.length][2];

        for (int col = 0; col < array[0].length; col++) {

            int counter = 0;

            for (int row = 0; row < array.length; row++) {

                if (!(((Double)array[row][col]).isNaN())) {
                    counter += 1;
                }
            }
            double[] columnNoNaNArray = new double[counter];

            int index = 0;
            for (int row = 0; row < array.length; row++) {

                if (!(((Double)array[row][col]).isNaN())) {
                    columnNoNaNArray[index] = array[row][col];
                    index += 1;
                }
            }
            int[] indices = new int[counter];
            HeapSort.sort(columnNoNaNArray, indices);
            double hpdBinInterval[] = getHPDInterval(0.95, columnNoNaNArray, indices);

            returnArray[col][0] = hpdBinInterval[0];
            returnArray[col][1] = hpdBinInterval[1];

        }

    return returnArray;
    }

}
