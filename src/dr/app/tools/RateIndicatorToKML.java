package dr.app.tools;

import dr.evolution.io.Importer;
import dr.util.HeapSort;
import dr.stats.DiscreteStatistics;

import java.io.*;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: phil
 * Date: Feb 14, 2008
 * Time: 8:13:35 AM
 * To change this template use File | Settings | File Templates.
 *
 *  RateIndicatorLog locationDegrees burnin poissonMean poissonOffset cutoff
 *
 */
public class RateIndicatorToKML {
    public static void main(String[] args) throws java.io.IOException, Importer.ImportException {

        String RateIndicatorLog = args[0];
        String rateIndicatorString	= "rateIndicator";
        String actualRateString = "productStatistic";
        String relativeRateString = "relativeRate";
        String geoSiteModelString = "geoSiteModel";
        boolean nonreversible = false;
        boolean truncated = true; // a nontruncated is used for nonreversible matrices

        String StatesAndCoordinates = args[1];

        int burnin = Integer.parseInt(args[2]);
        System.out.println("Ignoring first " + burnin + " trees as burnin.");

        double meanPoissonPrior = Double.parseDouble(args[3]);
        int offsetPoissonPrior = Integer.parseInt(args[4]);
        if (!(offsetPoissonPrior > 0)) {
            truncated = false;
            System.out.println("Poisson prior not truncated, offset = " + offsetPoissonPrior);
        }

        String cutoff = args[5];
        double cutoffValue;
        boolean BFcutoff = false;
        if (cutoff.contains("bf")) {
           BFcutoff = true;
           cutoffValue = Double.parseDouble(cutoff.replace("bf",""));
        }  else {
           BFcutoff = false;
           cutoffValue = Double.parseDouble(cutoff);
        }

        String startColorHighest = "E6F0FA"; //startpoint is the youngest node
        String endColorHighest = "9314FF"; //

        String startColorMRC = "00FF00"; //startpoint is the youngest node
        String endColorMRC = "00F1D6"; //

        //count the number of states in the RateIndicatorLog file
        BufferedReader reader1 = new BufferedReader(new FileReader(RateIndicatorLog));
        String current1 = reader1.readLine();
        int counter1 = 0;
        while (current1 != null && !reader1.equals("")) {
            while (current1.startsWith("#")) {
                current1 = reader1.readLine();
            }
            current1 = reader1.readLine();
            counter1++;
        }
        System.out.println("number of generations: "+(counter1 - 1));

        //find the first rateIndicator in the RateIndicatorLog file
        int firstRateIndicator = getFirstEntryOf(RateIndicatorLog, rateIndicatorString);
        System.out.println("first rateIndicator is at column "+firstRateIndicator);

        //boolean to see it the rateLOg contains productStatitistics, if not, we make 'em ourselves
        boolean containsActualRates = hasEntryOf(RateIndicatorLog, actualRateString);
        if (!containsActualRates){
            actualRateString = relativeRateString;
        }

        int firstActualRate = getFirstEntryOf(RateIndicatorLog, actualRateString);
        System.out.println("first actual rate (productStatistic) is at column "+firstActualRate);

        int geoSiteModelMuPosition = getFirstEntryOf(RateIndicatorLog, geoSiteModelString);
        System.out.println("geoSiteModelMu is at column "+geoSiteModelMuPosition);

        //count the rateIndicators in the RateIndicatorLog file
        int numberOfRateIndicators = getNumberOfEntries(RateIndicatorLog, firstRateIndicator, rateIndicatorString);
        int numberOfActualRates = numberOfRateIndicators;

        // check whether the number of numberOfrateIndicators matches the number of states
        BufferedReader readerStates = new BufferedReader(new FileReader(StatesAndCoordinates));
        String currentStates = readerStates.readLine();
        int statesCounter = 0;
        while (currentStates != null && !readerStates.equals("")) {
            currentStates = readerStates.readLine();
            statesCounter++;
        }
        System.out.println("number of states in "+StatesAndCoordinates+" = "+statesCounter);
        if (numberOfRateIndicators != ((statesCounter*(statesCounter-1.0))/2.0)) {
              if (numberOfRateIndicators == (statesCounter*(statesCounter-1.0))) {
                System.out.println("K*(K-1) rateIndicators, with K = "+ statesCounter+". So, nonreversible matrix!");
                nonreversible = true;
              } else {
              System.err.println("the number of rateIndicators ("+numberOfRateIndicators+") does not match (K*(K-1)/2) states ("+((statesCounter*(statesCounter-1.0))/2.0)+"; with K = "+statesCounter+").");
              }
        }

        // so now we know the dimension of the rateIndicator array
        double[][] rateIndicators = new double[((counter1 - 1)-burnin)][numberOfRateIndicators];

        // now read rateIndicators in rateIndicator array
        BufferedReader indicatorReader = new BufferedReader(new FileReader(RateIndicatorLog));
        String rateIndicatorCurrent = indicatorReader.readLine();
        while (rateIndicatorCurrent.startsWith("#")) {
            rateIndicatorCurrent = indicatorReader.readLine();
        }

        // skip the headers in the rateIndicator file
        while (rateIndicatorCurrent.startsWith("state")) {
            rateIndicatorCurrent = indicatorReader.readLine();
        }

        double rateIndicator;
        int linesRead = 0;
        //skip burnin
        while (linesRead < burnin){
            rateIndicatorCurrent = indicatorReader.readLine();
            linesRead ++;
        }

        int rowCounter = 0;

        while (rateIndicatorCurrent!= null && !indicatorReader.equals("")) {

            int columnCounter = 0;

            int startCounter = 1;
            int rateIndicatorCounter = 0;

            StringTokenizer tokens = new StringTokenizer(rateIndicatorCurrent);
            rateIndicator = Double.parseDouble(tokens.nextToken());

            //skip until we encounter rateIndicators
            while (startCounter < firstRateIndicator) {
                rateIndicator = Double.parseDouble(tokens.nextToken());
                startCounter ++;

            }

            // read all rateIndicators
            while (rateIndicatorCounter < numberOfRateIndicators) {
                rateIndicators[rowCounter][columnCounter] = rateIndicator;
                columnCounter ++;
                rateIndicatorCounter ++;
                rateIndicator = Double.parseDouble(tokens.nextToken());

            }

            rowCounter ++;
            rateIndicatorCurrent = indicatorReader.readLine();

        }  // rateIndicators are now in the rateIndicator array
        //print2DArray(rateIndicators, "test.txt");

        double[] expectedRateIndicators = meanCol(rateIndicators);

        // read in locations and degrees
        String[] locations = new String[statesCounter];
        double[] longitudes = new double[statesCounter];
        double[] latitudes = new double[statesCounter];
        BufferedReader locationsReader = new BufferedReader(new FileReader(StatesAndCoordinates));
        String currentLocation = locationsReader.readLine();
        int locationCounter = 0;
        while (currentLocation != null && !locationsReader.equals("")) {
            StringTokenizer tokens = new StringTokenizer(currentLocation);
            locations[locationCounter] = tokens.nextToken();
            System.out.print(locations[locationCounter]+"\t");
            latitudes[locationCounter] = Double.parseDouble(tokens.nextToken());
            System.out.print(latitudes[locationCounter]+"\t");
            longitudes[locationCounter] = Double.parseDouble(tokens.nextToken());
            System.out.print(longitudes[locationCounter]+"\r");
            currentLocation = locationsReader.readLine();
            locationCounter++;

        }

        // an array to get the actual rates and to get the geoSiteModel.mu's
        double[][] actualRates = new double[((counter1 - 1)-burnin)][numberOfActualRates];
        double[] geoSiteModelMus = new double[((counter1 - 1)-burnin)];

        // now lets get those rates and mu's
        BufferedReader dispersalReader = new BufferedReader(new FileReader(RateIndicatorLog));
        String dispersalCurrent = dispersalReader.readLine();
        while (dispersalCurrent.startsWith("#")) {
            dispersalCurrent = dispersalReader.readLine();
        }

        // skip the headers in the rateIndicator file
        while (dispersalCurrent.startsWith("state")) {
            dispersalCurrent = dispersalReader.readLine();
        }
        int linesRead2 = 0;
        //skip burnin
        while (linesRead2 < burnin){
            dispersalCurrent = dispersalReader.readLine();
            linesRead2 ++;
        }
        int rowCounter2 = 0;
        double actualRate;
        double geoSiteModelMu;
        while (dispersalCurrent!= null && !dispersalReader.equals("")) {

            int columnCounter = 0;

            int startRateCounter = 1;
            int actualRateCounter = 0;

            StringTokenizer tokens = new StringTokenizer(dispersalCurrent);
            actualRate = Double.parseDouble(tokens.nextToken());

            //skip until we encounter rateIndicators
            while (startRateCounter < firstActualRate) {
                actualRate = Double.parseDouble(tokens.nextToken());
                startRateCounter ++;
            }

            // compile expected distance per time unit
            while (actualRateCounter < numberOfRateIndicators) {
                actualRates[rowCounter2][columnCounter] = actualRate;
                columnCounter ++;
                actualRateCounter ++;
                actualRate = Double.parseDouble(tokens.nextToken());
            }

            StringTokenizer tokens2 = new StringTokenizer(dispersalCurrent);
            geoSiteModelMu = Double.parseDouble(tokens2.nextToken());
            int geoSiteModelMuCounter = 1;

            while (geoSiteModelMuCounter < geoSiteModelMuPosition) {
                geoSiteModelMu = Double.parseDouble(tokens2.nextToken());
                geoSiteModelMuCounter ++;
            }
            geoSiteModelMus[rowCounter2] = geoSiteModelMu;
            //System.out.println(geoSiteModelMus[rowCounter2]);

            rowCounter2 ++;
            dispersalCurrent = dispersalReader.readLine();

        }


       //if we hadn't read in actual rates (productStatistics =  relative rate * rateIndicator), we make those here
        if (!containsActualRates){
            for (int a = 0; a < actualRates.length; a++){
                for (int b = 0; b < actualRates[0].length; b++){
                    actualRates[a][b] = actualRates[a][b]*rateIndicators[a][b];
                }
            }
        }


        // array for expected distances per unit time
        double[] expectedDistancesPerUnitTime = new double[(counter1 - 1)-burnin];
        double[] distances = new double[numberOfActualRates];
        if (!nonreversible) {
            int pairwiseComparisons = 0;
            int pairwiseCounter1 = 0;
            for (int c = 0; c < statesCounter; c++) {
                pairwiseCounter1 ++;
                for (int d = pairwiseCounter1; d < statesCounter; d++) {

                    double lat1 = latitudes[c];
                    double lat2 = latitudes[d];
                    double lon1 = longitudes[c];
                    double lon2 = longitudes[d];
                    distances[pairwiseComparisons] = (3958*Math.PI*Math.sqrt((lat2-lat1)*(lat2-lat1)+Math.cos(lat2/57.29578)*Math.cos(lat1/57.29578)*(lon2-lon1)*(lon2-lon1))/180);  //3958 is the intermediate of 3963 miles at the equator versus 3950 miles at the poles
                    pairwiseComparisons ++;

                }
            }
        } else {
            int distanceCounter = 0;
            for (int a = 0; a < statesCounter; a++) {
                for (int b = 0; b < statesCounter; b++) {
                    double lat1 = latitudes[a];
                    double lat2 = latitudes[b];
                    double lon1 = longitudes[a];
                    double lon2 = longitudes[b];
                    double distance = (3958*Math.PI*Math.sqrt((lat2-lat1)*(lat2-lat1)+Math.cos(lat2/57.29578)*Math.cos(lat1/57.29578)*(lon2-lon1)*(lon2-lon1))/180);
                    if (distance > 0) {
                        distances[distanceCounter] = distance;
                        distanceCounter ++;
                    }
                }
            }

        }

        // now let's get those expected distances per unit time
        for (int a = 0; a < expectedDistancesPerUnitTime.length; a++) {
            double sumOfRateDistanceFrequency = 0.0;
            for (int b = 0; b < numberOfActualRates; b ++) {
                sumOfRateDistanceFrequency += actualRates[a][b]*(1.0/statesCounter)*distances[b];
            }
            expectedDistancesPerUnitTime[a] = geoSiteModelMus[a]*sumOfRateDistanceFrequency;
            //System.out.println(expectedDistancesPerUnitTime[a]);

        }

        // mean and HPDs for the expected distance per unit time
        double meanDistancePerUnitTime = DiscreteStatistics.mean(expectedDistancesPerUnitTime);
        double medianDistancePerUnitTime = DiscreteStatistics.quantile(0.5, expectedDistancesPerUnitTime);
        int[] indicesHPD = new int[expectedDistancesPerUnitTime.length];
        HeapSort.sort(expectedDistancesPerUnitTime, indicesHPD);
        double hpdInterval[] = getHPDInterval(0.95, expectedDistancesPerUnitTime, indicesHPD);
        System.out.println("\rmean distance per time unit\tmedian\tHPDlow\tHPDup");
        System.out.println(meanDistancePerUnitTime+"\t"+medianDistancePerUnitTime+"\t"+hpdInterval[0]+"\t"+hpdInterval[1]+"\r\r");

        //
        BufferedWriter dispersalRate = new BufferedWriter(new FileWriter("dispersalRates.txt"));
        StringBuffer dispersalRateBuffer = new StringBuffer();
        dispersalRateBuffer.append("state\tdispersalRate(miles/year)\tdispersalRate(km/year)\r");
        for (int s = 0; s < expectedDistancesPerUnitTime.length; s++){
            dispersalRateBuffer.append(s+"\t"+expectedDistancesPerUnitTime[s]+"\t"+(expectedDistancesPerUnitTime[s]*1.609344)+"\r");
        }
        dispersalRate.write(dispersalRateBuffer.toString());
        dispersalRate.close();

        //write locations to KML
        StringBuffer locationBuffer = new StringBuffer();
        for (int x = 0; x < statesCounter; x++) {
            writeLocationToKML(locations[x], longitudes[x], latitudes[x], locationBuffer);
        }


        // make arrays with start and end locations/coordinates
        String[] startLocations = new String [numberOfRateIndicators];
        String[] endLocations = new String [numberOfRateIndicators];
        double[] startLongitude = new double [numberOfRateIndicators];
        double[] endLongitude = new double [numberOfRateIndicators];
        double[] startLatitude = new double [numberOfRateIndicators];
        double[] endLatitude = new double [numberOfRateIndicators];

        /** this was for the old rate ordering in the nonreversible, obsolete now
        if (!nonreversible) {
            int elementCounter = 0;
            int secondCounter = 0;
            for (int i = 0; i < (statesCounter - 1); i++) {

                secondCounter ++;

                for (int j = secondCounter; j < statesCounter; j++) {

                    startLocations[elementCounter] = locations[i];
                    startLongitude[elementCounter] = longitudes[i];
                    startLatitude[elementCounter] = latitudes[i];
                    endLocations[elementCounter] = locations[j];
                    endLongitude[elementCounter] = longitudes[j];
                    endLatitude[elementCounter] = latitudes[j];

                    elementCounter ++;
                }
            }
        } else {
            int elementCounter = 0;
            for (int i = 0; i < statesCounter; i++) {

                for (int j = 0; j < statesCounter; j++) {

                    if (i != j) {

                        startLocations[elementCounter] = locations[i];
                        startLongitude[elementCounter] = longitudes[i];
                        startLatitude[elementCounter] = latitudes[i];
                        endLocations[elementCounter] = locations[j];
                        endLongitude[elementCounter] = longitudes[j];
                        endLatitude[elementCounter] = latitudes[j];

                        elementCounter ++;

                    }
                }
            }

        }
         end of older code **/

        //begin of new code
        int elementCounter = 0;
        int secondCounter = 0;
        for (int i = 0; i < (statesCounter - 1); i++) {

            secondCounter ++;

                for (int j = secondCounter; j < statesCounter; j++) {

                    startLocations[elementCounter] = locations[i];
                    startLongitude[elementCounter] = longitudes[i];
                    startLatitude[elementCounter] = latitudes[i];
                    endLocations[elementCounter] = locations[j];
                    endLongitude[elementCounter] = longitudes[j];
                    endLatitude[elementCounter] = latitudes[j];

                    elementCounter ++;
                }
        }
        // for nonreversible models, we keep on filling the arrays
        if (nonreversible) {
               for (int k = 0; k < elementCounter; k++) {

                   startLocations[elementCounter + k] = endLocations[k];
                   startLongitude[elementCounter + k] = endLongitude[k];
                   startLatitude[elementCounter + k] = endLatitude[k];
                   endLocations[elementCounter + k] = startLocations[k];
                   endLongitude[elementCounter + k] = startLongitude[k];
                   endLatitude[elementCounter + k] = startLatitude[k];
               }
        }
        //end of new code **/

        //print actualrates to csv file for nonreversible rates
        double[] meanRates = meanCol(actualRates);
        if (nonreversible) {
            printRateMatrix(meanRates, locations, startLocations, endLocations, "nonreversible.csv");
        }

        //sort expected rateIndicator
        int[] indices = new int[numberOfRateIndicators];
        HeapSort.sort(expectedRateIndicators, indices);

        //print expected rateIndicators
        for (int a = 0; a < numberOfRateIndicators; a++) {

            System.out.println(
                expectedRateIndicators[indices[a]]+"\tmeanRate="+meanRates[indices[a]]+": between "+startLocations[indices[a]]+" (long: "+startLongitude[indices[a]]+"; lat: "+ startLatitude[indices[a]]+")" +
                " and "+ endLocations[indices[a]]+" (long: "+endLongitude[indices[a]]+"; lat: "+ endLatitude[indices[a]]+")"
            );
        }

        double branchWidthConstant = 1.0;
        double branchWidthMultiplier = 6.0;

        StringBuffer rateBuffer = new StringBuffer();
        StringBuffer rateAltitudeBuffer = new StringBuffer();
        StringBuffer styleBuffer = new StringBuffer();

        System.out.println("\r\r");

        // write out the rateIndicators that are above the cutoff (real number or based on Bayes factor). The first time we write out to the screen and log the max and min rateIndicator
        double[] minAndMaxRateIndicator = new double[2];
        minAndMaxRateIndicator[0] = Double.MAX_VALUE;
        minAndMaxRateIndicator[1] = 0.0;

        if (BFcutoff) {
            System.out.println("Bayes factor cutoff = "+getBayesFactorCutOff(cutoffValue, meanPoissonPrior, offsetPoissonPrior, statesCounter, nonreversible));
            System.out.println(cutoffValue+"\t"+meanPoissonPrior+"\t"+statesCounter);

        }

        for (int o = 0; o < numberOfRateIndicators; o++){

            double indicator = expectedRateIndicators[indices[(numberOfRateIndicators - o - 1)]];

            if (!BFcutoff) {

                if (indicator > cutoffValue) {

                    System.out.println(
                        "I="+indicator+"\tBF="+getBayesFactor(indicator, meanPoissonPrior, statesCounter, offsetPoissonPrior, nonreversible)+": between "+startLocations[indices[(numberOfRateIndicators - o - 1)]]+" (long: "+startLongitude[indices[(numberOfRateIndicators - o - 1)]]+"; lat: "+ startLatitude[indices[(numberOfRateIndicators - o - 1)]]+")" +
                        " and "+ endLocations[indices[(numberOfRateIndicators - o - 1)]]+" (long: "+endLongitude[indices[(numberOfRateIndicators - o - 1)]]+"; lat: "+ endLatitude[indices[(numberOfRateIndicators - o - 1)]]+")"
                    );

                    if (indicator > minAndMaxRateIndicator[1]) {
                        minAndMaxRateIndicator[1] = indicator;
                    } else if (indicator < minAndMaxRateIndicator[0]) {
                        minAndMaxRateIndicator[0] = indicator;
                    }


                }

            }  else {

                if (indicator > getBayesFactorCutOff(cutoffValue, meanPoissonPrior, offsetPoissonPrior, statesCounter, nonreversible)) {

                    System.out.println(
                        "I="+indicator+"\tBF="+getBayesFactor(indicator, meanPoissonPrior, statesCounter, offsetPoissonPrior, nonreversible)+": between "+startLocations[indices[(numberOfRateIndicators - o - 1)]]+" (long: "+startLongitude[indices[(numberOfRateIndicators - o - 1)]]+"; lat: "+ startLatitude[indices[(numberOfRateIndicators - o - 1)]]+")" +
                        " and "+ endLocations[indices[(numberOfRateIndicators - o - 1)]]+" (long: "+endLongitude[indices[(numberOfRateIndicators - o - 1)]]+"; lat: "+ endLatitude[indices[(numberOfRateIndicators - o - 1)]]+")"
                    );

                    if (indicator > minAndMaxRateIndicator[1]) {
                        minAndMaxRateIndicator[1] = indicator;
                    } else if (indicator < minAndMaxRateIndicator[0]) {
                        minAndMaxRateIndicator[0] = indicator;
                    }


                }
            }
        }

        double altitudeFactor = 200;
        double divider = 100;

        //outputfile for GenGis
        BufferedWriter out2 = new BufferedWriter(new FileWriter("connections.csv"));
        out2.append("Node1,Node2,Weight");
        out2.newLine();

        int indicatorCounter2 = 0;
        for (int p = 0; p < numberOfRateIndicators; p++){

            double indicator = expectedRateIndicators[indices[(numberOfRateIndicators - p - 1)]];

            if (!BFcutoff) {

                if (indicator > cutoffValue) {

                    //writeRatesToKML(indicator, startLongitude[indices[(numberOfRateIndicators - p - 1)]], startLatitude[indices[(numberOfRateIndicators - p - 1)]], endLongitude[indices[(numberOfRateIndicators - p - 1)]], endLatitude[indices[(numberOfRateIndicators - p - 1)]], (indicatorCounter2 + 1), branchWidthConstant, branchWidthMultiplier, minAndMaxRateIndicator, rateBuffer, styleBuffer);
                    writeRatesWithAltitudeToKML(indicator, startLongitude[indices[(numberOfRateIndicators - p - 1)]], startLatitude[indices[(numberOfRateIndicators - p - 1)]], endLongitude[indices[(numberOfRateIndicators - p - 1)]], endLatitude[indices[(numberOfRateIndicators - p - 1)]], (indicatorCounter2 + 1), branchWidthConstant, branchWidthMultiplier, minAndMaxRateIndicator, altitudeFactor, divider, startColorHighest, endColorHighest, rateAltitudeBuffer, styleBuffer);
                    out2.append(startLocations[indices[(numberOfRateIndicators - p - 1)]]+","+endLocations[indices[(numberOfRateIndicators - p - 1)]]+","+indicator);
                    out2.newLine();
                    indicatorCounter2 ++;

                }

            }  else {

                if (indicator > getBayesFactorCutOff(cutoffValue, meanPoissonPrior, offsetPoissonPrior, statesCounter, nonreversible)) {

                    //writeRatesToKML(indicator, startLongitude[indices[(numberOfRateIndicators - p - 1)]], startLatitude[indices[(numberOfRateIndicators - p - 1)]], endLongitude[indices[(numberOfRateIndicators - p - 1)]], endLatitude[indices[(numberOfRateIndicators - p - 1)]], (indicatorCounter2 + 1), branchWidthConstant, branchWidthMultiplier, minAndMaxRateIndicator, rateBuffer, styleBuffer);
                    writeRatesWithAltitudeToKML(indicator, startLongitude[indices[(numberOfRateIndicators - p - 1)]], startLatitude[indices[(numberOfRateIndicators - p - 1)]], endLongitude[indices[(numberOfRateIndicators - p - 1)]], endLatitude[indices[(numberOfRateIndicators - p - 1)]], (indicatorCounter2 + 1), branchWidthConstant, branchWidthMultiplier, minAndMaxRateIndicator, altitudeFactor, divider, startColorHighest, endColorHighest, rateAltitudeBuffer, styleBuffer);
                    out2.append(startLocations[indices[(numberOfRateIndicators - p - 1)]]+","+endLocations[indices[(numberOfRateIndicators - p - 1)]]+","+indicator);
                    out2.newLine();
                    indicatorCounter2 ++;

                }
            }
        }

        // lets now find the maximum rate credibility configuration.
        double maxRateCredibility = Double.MAX_VALUE;
        int indexOfmaxRateCredibility = -1;
        for (int q = 0; q < rateIndicators.length; q++) {

            double logScore = 0.0;

            for (int r = 0; r < rateIndicators[0].length; r++) {

                if (rateIndicators[q][r] > 0) {
                    logScore += Math.log(expectedRateIndicators[r]);
                }


            }

            if (Math.abs(logScore) < maxRateCredibility) {

                indexOfmaxRateCredibility = q;
                maxRateCredibility = Math.abs(logScore);

            }
        }

        System.out.println("max rate credibity = -"+maxRateCredibility+", at index = "+indexOfmaxRateCredibility);


        //now that we have found the maximum rate credibility configuration, plot it as kml
        StringBuffer maxRateBuffer = new StringBuffer();
        double[] minAndMaxRate = new double[2];
        minAndMaxRate[0] = Double.MAX_VALUE;
        minAndMaxRate[1] = 0.0;
        for (int s = 0; s < rateIndicators[indexOfmaxRateCredibility].length; s++) {

            if (rateIndicators[indexOfmaxRateCredibility][s] > 0) {

                if (expectedRateIndicators[s] < minAndMaxRate[0]) {

                    minAndMaxRate[0] = expectedRateIndicators[s];

                } else if (expectedRateIndicators[s] > minAndMaxRate[1]) {

                    minAndMaxRate[1] = expectedRateIndicators[s];

                }

                System.out.println(
                    expectedRateIndicators[s]+": between "+startLocations[s]+" (long: "+startLongitude[s]+"; lat: "+ startLatitude[s]+")" +
                    " and "+ endLocations[s]+" (long: "+endLongitude[s]+"; lat: "+ endLatitude[s]+")"
                );


            }
        }
        for (int t = 0; t < rateIndicators[indexOfmaxRateCredibility].length; t++) {

            if (rateIndicators[indexOfmaxRateCredibility][t] > 0) {

                //we give it the number (numberOfRateIndicators + t + 1) to have different numbers than the highest rateIndicators
                writeRatesToKML(expectedRateIndicators[t], startLongitude[t], startLatitude[t], endLongitude[t], endLatitude[t], (numberOfRateIndicators + t + 1), branchWidthConstant, branchWidthMultiplier, minAndMaxRate, startColorMRC, endColorMRC,  maxRateBuffer, styleBuffer);

            }
        }



        BufferedWriter out1 = new BufferedWriter(new FileWriter("rates.kml"));
        StringBuffer buffer = new StringBuffer();
        compileKMLBuffer(locationBuffer, rateBuffer, rateAltitudeBuffer, maxRateBuffer, styleBuffer, buffer);
        out1.write(buffer.toString());
        out1.close();
        out2.close();
    }



    private static int getFirstEntryOf(String file, String entryString)    {

        int firstRateIndicator = 1;

        try {
            BufferedReader reader2 = new BufferedReader(new FileReader(file));
            String current2 = reader2.readLine();

            //skip comment lines
            while (current2.startsWith("#")) {
                current2 = reader2.readLine();
            }

            String parameter1 = null;
            StringTokenizer tokens1 = new StringTokenizer(current2);
            parameter1 = tokens1.nextToken();
            while(!parameter1.contains(entryString)) {
                parameter1 = tokens1.nextToken();
                firstRateIndicator ++;
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return firstRateIndicator;

    }

    private static boolean hasEntryOf(String file, String entryString)    {

        boolean hasEntry = false;

        try {
            BufferedReader reader2 = new BufferedReader(new FileReader(file));
            String current2 = reader2.readLine();

            //skip comment lines
            while (current2.startsWith("#")) {
                current2 = reader2.readLine();
            }

            String parameter1 = null;
            StringTokenizer tokens1 = new StringTokenizer(current2);

            while(tokens1.hasMoreTokens()) {
                parameter1 = tokens1.nextToken();
                if (parameter1.contains(entryString)) {
                    hasEntry = true;    
                }
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return hasEntry;

    }

    private static int getNumberOfEntries(String file, int firstRateIndicator, String rateIndicatorString)    {

        int numberOfRateIndicators = 0;

        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String current3 = reader.readLine();

            //skip comment lines
            while (current3.startsWith("#")) {
                current3 = reader.readLine();
            }

            String rateIndicator = null;
            int startCounter = 1;
            StringTokenizer tokens = new StringTokenizer(current3);
            rateIndicator = tokens.nextToken();
            //find first rateIndicator..
            while (startCounter < firstRateIndicator) {
                  rateIndicator = tokens.nextToken();
                  startCounter ++;
            }
            // and continue counting from thereon
            while(rateIndicator.contains(rateIndicatorString)) {
                  rateIndicator = tokens.nextToken();
                  numberOfRateIndicators ++;
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return numberOfRateIndicators;

    }

    private static double[] meanCol(double[][] x)    {
         double[] returnArray = new double[x[0].length];

         for (int i = 0; i < x[0].length; i++) {

             double m = 0.0;
             int len = 0;

             for (int j = 0; j < x.length; j++) {

                m += x[j][i];
                len += 1;

             }

             returnArray[i] = m / (double) len;

         }
         return returnArray;
     }

    private static double getBayesFactorCutOff(double bayesFactor, double meanPoissonPrior, int offset, int numberOfLocations, boolean nonreversible) {

        double bayesFactorCutoff = 0;
        double posteriorOdds = 0;
        int numberOfRatesMultiplier = 2;
        if (nonreversible) numberOfRatesMultiplier = 1;

        double priorProbability = (meanPoissonPrior + offset)/((numberOfLocations*(numberOfLocations-1))/numberOfRatesMultiplier);
        double priorOdds = priorProbability/(1.0 - priorProbability);
        posteriorOdds = priorOdds*bayesFactor;
        bayesFactorCutoff = posteriorOdds/(1.0 + posteriorOdds);

        return bayesFactorCutoff;

    }

    private static double getBayesFactor(double meanIndicator, double meanPoissonPrior, int numberOfLocations, int offset, boolean nonreversible) {

        double bayesFactor = 0;
        double priorProbabilityDenominator = 0;
        int numberOfRatesMultiplier = 2;
        priorProbabilityDenominator = meanPoissonPrior + offset;
        if (nonreversible) numberOfRatesMultiplier = 1;

        double priorProbability = priorProbabilityDenominator/((numberOfLocations*(numberOfLocations-1))/numberOfRatesMultiplier);

        double priorOdds = priorProbability/(1.0 - priorProbability);
        double posteriorProbability = meanIndicator;
        double posteriorOdds =  posteriorProbability/(1 - posteriorProbability);
        bayesFactor = posteriorOdds/priorOdds;

        return bayesFactor;

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

    private static void writeRatesToKML(double rateIndicator, double startLongitude, double startLatitude, double endLongitude, double endLatitude, int number, double branchWidthConstant, double branchWidthMultiplier, double[] minAndMaxRateIndicator, String startColor, String endColor, StringBuffer rateBuffer, StringBuffer styleBuffer) {

        rateBuffer.append("\t\t<Placemark>\r");
            rateBuffer.append("\t\t\t<name>rate"+number+"</name>\r");
            rateBuffer.append("\t\t\t<visibility>0</visibility>\r");
            rateBuffer.append("\t\t\t<styleUrl>#rate"+number+"_style</styleUrl>\r");
            rateBuffer.append("\t\t\t<LineString>\r");

                rateBuffer.append("\t\t\t\t<altitudeMode>clambToGround</altitudeMode>\r");
                rateBuffer.append("\t\t\t\t<tessellate>1</tessellate>\r");

                rateBuffer.append("\t\t\t\t<coordinates>\r");
                    rateBuffer.append("\t\t\t\t\t"+startLongitude+","+startLatitude+",0\r");
                    rateBuffer.append("\t\t\t\t\t"+endLongitude+","+endLatitude+",0\r");
                rateBuffer.append("\t\t\t\t</coordinates>\r");

            rateBuffer.append("\t\t\t</LineString>\r");
        rateBuffer.append("\t\t\t</Placemark>\r");

        styleBuffer.append("\t<Style id=\"rate"+number+"_style\">\r");
        styleBuffer.append("\t\t<LineStyle>\r");
        styleBuffer.append("\t\t\t<width>"+(branchWidthConstant+branchWidthMultiplier*rateIndicator)+"</width>\r");
        styleBuffer.append("\t\t\t<color>"+"FF"+getKMLColor(rateIndicator,minAndMaxRateIndicator, startColor, endColor)+"</color>\r");
        styleBuffer.append("\t\t</LineStyle>\r");
        styleBuffer.append("\t</Style>\r");

    }

       private static void writeRatesWithAltitudeToKML(double rateIndicator, double startLongitude, double startLatitude, double endLongitude, double endLatitude, int number, double branchWidthConstant, double branchWidthMultiplier, double[] minAndMaxRateIndicator, double altitudeFactor, double divider, String startColor, String endColor, StringBuffer rateAltitudeBuffer, StringBuffer styleBuffer) {

        double distance = (3958*Math.PI*Math.sqrt((endLatitude-startLatitude)*(endLatitude-startLatitude)+Math.cos(endLatitude/57.29578)*Math.cos(startLatitude/57.29578)*(endLongitude-startLongitude)*(endLongitude-startLongitude))/180);

        double maxAltitude = distance*altitudeFactor;
        double latitudeDifference = endLatitude - startLatitude;
        double longitudeDifference = endLongitude - startLongitude;


           boolean longitudeBreak = false; //if we go through the 180
           if (endLongitude*startLongitude < 0) {

               double trialDistance = 0;

               if (endLongitude < 0) {
                   trialDistance += (endLongitude + 180);
                   trialDistance += (180 - startLongitude);
                   //System.out.println(parentLongitude+"\t"+longitude+"\t"+trialDistance+"\t"+longitudeDifference);
               } else {
                   trialDistance += (startLongitude + 180);
                   trialDistance += (180 - endLongitude);
                   //System.out.println(parentLongitude+"\t"+longitude+"\t"+trialDistance+"\t"+longitudeDifference);
                }

               if (trialDistance < Math.abs(longitudeDifference)) {
                   longitudeDifference = trialDistance;
                   longitudeBreak = true;
                   //System.out.println("BREAK!"+longitudeDifference);
               }
           }
           

           styleBuffer.append("\t<Style id=\"rate"+number+"_style\">\r");
               styleBuffer.append("\t\t<LineStyle>\r");
                   styleBuffer.append("\t\t\t<width>"+(branchWidthConstant+branchWidthMultiplier*rateIndicator)+"</width>\r");
                   //styleBuffer.append("\t\t\t<color>"+"FF"+getWhiteToRedColor(rateIndicator,minAndMaxRateIndicator)+"</color>\r");
                   styleBuffer.append("\t\t\t<color>"+"FF"+getKMLColor(rateIndicator,minAndMaxRateIndicator, startColor, endColor)+"</color>\r");
           styleBuffer.append("\t\t</LineStyle>\r");
           styleBuffer.append("\t</Style>\r");

           double currentLongitude1 = 0;  //we need this if we go through the 180
           double currentLongitude2 = 0;  //we need this if we go through the 180


        for (int a = 0; a < divider; a ++) {


                rateAltitudeBuffer.append("\t\t<Placemark>\r");
                    rateAltitudeBuffer.append("\t\t\t<name>rate"+number+"_part"+(a+1)+"</name>\r");
                    rateAltitudeBuffer.append("\t\t\t<styleUrl>#rate"+number+"_style</styleUrl>\r");
                    rateAltitudeBuffer.append("\t\t\t<LineString>\r");
                        rateAltitudeBuffer.append("\t\t\t\t<altitudeMode>relativeToGround</altitudeMode>\r");
                        rateAltitudeBuffer.append("\t\t\t\t<tessellate>1</tessellate>\r");
                        rateAltitudeBuffer.append("\t\t\t\t<coordinates>\r");

                        if (longitudeBreak) {

                            if (startLongitude > 0) {
                                currentLongitude1 = startLongitude+a*(longitudeDifference/divider);

                                if (currentLongitude1 < 180) {
                                   rateAltitudeBuffer.append("\t\t\t\t\t"+currentLongitude1+",");
                                   //System.out.println("1 currentLongitude1 < 180\t"+currentLongitude1+"\t"+longitude);

                                } else {
                                   rateAltitudeBuffer.append("\t\t\t\t\t"+(-180-(180-currentLongitude1))+",");
                                   //System.out.println("2 currentLongitude1 > 180\t"+currentLongitude1+"\t"+(-180-(180-currentLongitude1))+"\t"+longitude);
                                }
                            } else {
                                currentLongitude1 = startLongitude-a*(longitudeDifference/divider);

                                if (currentLongitude1 > (-180)) {
                                    rateAltitudeBuffer.append("\t\t\t\t\t"+currentLongitude1+",");
                                    //System.out.println("currentLongitude1 > -180\t"+currentLongitude1+"\t"+longitude);
                                 } else {
                                    rateAltitudeBuffer.append("\t\t\t\t\t"+(180+(currentLongitude1+180))+",");
                                    //System.out.println("currentLongitude1 > -180\t"+(180+(currentLongitude1+180))+"\t"+longitude);
                                 }
                            }

                        } else {
                            rateAltitudeBuffer.append("\t\t\t\t\t"+(startLongitude+a*(longitudeDifference/divider))+",");
                        }

                        rateAltitudeBuffer.append((startLatitude+a*(latitudeDifference/divider))+","+(maxAltitude*Math.sin(Math.acos(1 - a*(1.0/(divider/2.0)))))+"\r");

                        if (longitudeBreak) {

                            if (startLongitude > 0) {
                                currentLongitude2 = startLongitude+(a+1)*(longitudeDifference/divider);

                                if (currentLongitude2 < 180) {
                                    rateAltitudeBuffer.append("\t\t\t\t\t"+(currentLongitude2)+",");
                                } else {
                                    rateAltitudeBuffer.append("\t\t\t\t\t"+(-180-(180-currentLongitude2))+",");
                                }
                            } else {
                                currentLongitude2 = startLongitude-(a+1)*(longitudeDifference/divider);

                                if (currentLongitude2 > (-180)) {
                                    rateAltitudeBuffer.append("\t\t\t\t\t"+currentLongitude2+",");
                                } else {
                                    rateAltitudeBuffer.append("\t\t\t\t\t"+(180+(currentLongitude2+180))+",");
                                }
                            }

                        } else {
                            rateAltitudeBuffer.append("\t\t\t\t\t"+(startLongitude+(a+1)*(longitudeDifference/divider))+",");
                        }

                        rateAltitudeBuffer.append((startLatitude+(a+1)*(latitudeDifference/divider))+","+(maxAltitude*Math.sin(Math.acos(1 - (a+1)*(1.0/(divider/2.0)))))+"\r");
                        rateAltitudeBuffer.append("\t\t\t\t</coordinates>\r");
                    rateAltitudeBuffer.append("\t\t\t</LineString>\r");
                rateAltitudeBuffer.append("\t\t</Placemark>\r");
        }
    }

    private static void writeLocationToKML(String location, double longitude, double latitude, StringBuffer locationBuffer) {

        locationBuffer.append("\t\t\t<Placemark>\r");

            locationBuffer.append("\t\t\t\t<name>"+location+"</name>\r");

            locationBuffer.append("\t\t\t\t<Point>\r");
                locationBuffer.append("\t\t\t\t\t<altitudeMode>relativeToGround</altitudeMode>\r");
                locationBuffer.append("\t\t\t\t\t<coordinates>"+longitude+","+latitude+","+"0</coordinates>\r");
            locationBuffer.append("\t\t\t\t</Point>\r");

        locationBuffer.append("\t\t\t</Placemark>\r");

    }

     private static void compileKMLBuffer(StringBuffer locationBuffer, StringBuffer rateBuffer, StringBuffer rateAltitudeBuffer, StringBuffer maxRateBuffer, StringBuffer styleBuffer, StringBuffer  buffer) {

         buffer.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r");
         buffer.append("<kml xmlns=\"http://earth.google.com/kml/2.2\">\r");

         buffer.append("<Document>\r");

            buffer.append("\t<name>rates</name>\r");
            buffer.append(styleBuffer);

            buffer.append("\t<Folder>\r");
            buffer.append("\t\t<name>locations</name>\r");
            buffer.append(locationBuffer);
            buffer.append("\t</Folder>\r");
            /**
            buffer.append("\t<Folder>\r");
            buffer.append("\t\t<name>highest rateIndicators</name>\r");
            buffer.append("\t\t<description>highest rateIndicators based on expected rateIndicator value</description>\r");
            buffer.append(rateBuffer);
            buffer.append("\t</Folder>\r");  **/

             buffer.append("\t<Folder>\r");
             buffer.append("\t\t<name>highest rateIndicators</name>\r");
             buffer.append("\t\t<description>rates are shown for the highest expected rateIndicator value, either the ones above a Bayes factor cut-off or the highest 'x' rateIdicators</description>\r");
             buffer.append(rateAltitudeBuffer);
             buffer.append("\t</Folder>\r");

             buffer.append("\t<Folder>\r");
             buffer.append("\t\t<name>maximum rateIndicator credibility configuration</name>\r");
             buffer.append(maxRateBuffer);
             buffer.append("\t</Folder>\r");


         buffer.append("</Document>\r");
         buffer.append("</kml>\r");

     }

    private static String getBlueToRedColor(double value, double[] maxAndMin) {

        String color = "ffffffff";

        double colorUnits = 0;
        int colorInt = 0;

        colorUnits = ((maxAndMin[1] - value)/(maxAndMin[1]-maxAndMin[0]))*9.0;

        colorInt = (int)(Math.round(10 - colorUnits));

        switch (colorInt) {
            case 1:  color="B36600"; break;
            case 2:  color="CC3300"; break;
            case 3:  color="B31919"; break;
            case 4:  color="990033"; break;
            case 5:  color="990040"; break;
            case 6:  color="990066"; break;
            case 7:  color="990099"; break;
            case 8:  color="9900CC"; break;
            case 9:  color="6600E6"; break;
            case 10: color="0000FF"; break;
        }

        return color;

    }
    private static String getWhiteToRedColor(double value, double[] maxAndMin) {

        String color = "ffffffff";

        double colorUnits = 0;
        int colorInt = 0;

        colorUnits = ((maxAndMin[1] - value)/(maxAndMin[1]-maxAndMin[0]))*9.0;

        colorInt = (int)(Math.round(10 - colorUnits));

        switch (colorInt) {
            case 1:  color="FFFFFF"; break;
            case 2:  color="E5E5FF"; break;
            case 3:  color="CCCCFF"; break;
            case 4:  color="B2B2FF"; break;
            case 5:  color="9999FF"; break;
            case 6:  color="7F7FFF"; break;
            case 7:  color="6666FF"; break;
            case 8:  color="4C4CFF"; break;
            case 9:  color="3333FF"; break;
            case 10: color="1919FF"; break;
        }

        return color;


    }

    private static String getKMLColor(double value, double[] minMax, String startColor, String endColor) {

        startColor = startColor.toLowerCase();
        String startBlue = startColor.substring(0,2);
        String startGreen = startColor.substring(2,4);
        String startRed = startColor.substring(4,6);

        endColor =  endColor.toLowerCase();
        String endBlue = endColor.substring(0,2);
        String endGreen = endColor.substring(2,4);
        String endRed = endColor.substring(4,6);

        double proportion = (value - minMax[0])/(minMax[1] - minMax[0]);

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

        String color = blue+green+red;

        return color;
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

    private static void printRateMatrix(double[] rateArray, String[] locations, String[] startLocations, String[] endLocations, String name) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

            outFile.print("location");
            for (int o = 0; o < locations.length; o++) {
                outFile.print(","+locations[o]);
            }
            outFile.print("\r");

            for(int a = 0; a < locations.length; a++) {

                outFile.print(locations[a]);

                for(int b = 0; b < locations.length; b++) {

                    if (a == b) {

                        outFile.print(",-");

                    }  else {

                           for (int c = 0; c < rateArray.length; c ++) {

                            if ((startLocations[c].equals(locations[a]))&&(endLocations[c].equals(locations[b]))) {
                                outFile.print(","+rateArray[c]);
                            }

                           }

                    }

                }
                outFile.print("\r");
            }
/**            int location = 0;
            for (int p = 0; p < rateArray.length; p++) {

                if (p == 0) {
                    outFile.print(locations[location]);
                    location++;
                    outFile.print("\t-");
                }
                outFile.print("\t"+rateArray[p]);

                if (!((p+1)%(locations.length) > 0)) {
                    outFile.print("\t-");
                }
                if (!((p+1)%(locations.length-1) > 0)) {
                    outFile.print("\r");
                    if (location < locations.length){
                        outFile.print(locations[location]);
                        location++;
                    }
                }

            }
**/
            outFile.close();

        } catch(IOException io) {
           System.err.print("Error writing to file: " + name);
        }
    }

}

