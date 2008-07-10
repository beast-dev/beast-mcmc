package dr.evomodel.tree;

import dr.app.beast.BeastVersion;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodelxml.LoggerParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.logging.Logger;

/**
 * Projects and logs the random grid of the skyride onto a fixed grid.
 *
 * @author Erik Bloomquist
 */

public class GMRFSkyrideFixedGridLogger extends MCLogger {

    public static final String SKYRIDE_FIXED_GRID = "logSkyrideFixedGrid";
    public static final String NUMBER_OF_INTERVALS = "numberOfIntervals";
    public static final String GRID_HEIGHT = "gridHeight";
    public static final String AUXILLIARY_INTERVALS = "auxiliaryIntervals";


    private GMRFSkyrideLikelihood gmrfLike;
    private double[] fixedPopSize;
    private double intervalNumber;
    private double gridHeight;
    private double distance;
    private int auxiliaryIntervals;


    public GMRFSkyrideFixedGridLogger(LogFormatter formatter, int logEvery,
                                      GMRFSkyrideLikelihood gmrfLike,
                                      int intervalNumber, double gridHeight,
                                      int auxiliaryIntervals) {
        super(formatter, logEvery, false);
        this.gmrfLike = gmrfLike;
        this.intervalNumber = intervalNumber;
        this.gridHeight = gridHeight;
        this.distance = gridHeight / intervalNumber;
        this.auxiliaryIntervals = auxiliaryIntervals;

        fixedPopSize = new double[intervalNumber];
    }


    private void calculateGrid() {
        //The first random interval is the coalescent interval closest to the sampling time
        //The last random interval corresponds to the root.

        double[] treeInterval = gmrfLike.getCopyOfCoalescentIntervals();
        double[] treePopSize = gmrfLike.getPopSizeParameter().getParameterValues();

        double[] randomInterval = new double[treeInterval.length + auxiliaryIntervals];
        double[] randomPopSize = new double[randomInterval.length];

        for (int i = 1; i < treeInterval.length; i++)
            treeInterval[i] += treeInterval[i - 1];

        System.arraycopy(treeInterval, 0, randomInterval, 0, treeInterval.length);
        System.arraycopy(treePopSize, 0, randomInterval, 0, treePopSize.length);

        double length = (gridHeight - treeInterval[treeInterval.length - 1]) / auxiliaryIntervals;

        for (int i = treeInterval.length; i < randomInterval.length; i++) {

        }


        double intervalLength = gridHeight / intervalNumber;

        double startTime = 0;
        double endTime = 0;
        int interval = 0;
        for (int i = 0; i < fixedPopSize.length; i++) {
            fixedPopSize[i] = 0;
            if (interval < treeInterval.length) {

                endTime += intervalLength;
                double timeLeft = intervalLength;

                while (interval < treeInterval.length &&
                        treeInterval[interval] <= endTime) {
                    fixedPopSize[i] += (treeInterval[interval] - startTime - intervalLength + timeLeft) *
                            treePopSize[interval];
                    timeLeft = (intervalLength - treeInterval[interval] + startTime);
                    interval++;
                }

                if (interval < treeInterval.length) {
                    fixedPopSize[i] += timeLeft * treePopSize[interval];
                    fixedPopSize[i] /= intervalLength;
                } else {
                    fixedPopSize[i] /= (intervalLength - timeLeft);
                }

                startTime += intervalLength;

            } else {
                fixedPopSize[i] = -99;
            }
        }

    }

    private void getStringGrid(String[] values) {
        calculateGrid();

        int firstNA;
        for (firstNA = 0; firstNA < fixedPopSize.length; firstNA++) {
            if (fixedPopSize[firstNA] == -99) {
                break;
            }
        }

        for (int i = 0; i < firstNA; i++) {
            values[i + 1] = Double.toString(fixedPopSize[i]);
        }

        double sqrtPrecisionPerDistance = Math.sqrt((double) distance / gmrfLike.getPrecisionParameter().getParameterValue(0));

        double previous = fixedPopSize[firstNA - 1];

        for (int i = firstNA; i < fixedPopSize.length; i++) {
            double a = MathUtils.nextGaussian();
            double b = previous + a * sqrtPrecisionPerDistance;
            previous = b;
            values[i + 1] = Double.toString(b);
        }
        values[values.length - 1] = "" + -99;
    }

    public void log(int state) {
        if (logEvery <= 0 || ((state % logEvery) == 0)) {
            String[] values = new String[fixedPopSize.length + 1];

            values[0] = Integer.toString(state);

            getStringGrid(values);

            logValues(values);
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return "Projects the random grid of the skyride onto a fixed grid";
        }

        public Class getReturnType() {
            return GMRFSkyrideFixedGridLogger.class;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(NUMBER_OF_INTERVALS),
                AttributeRule.newDoubleRule(GRID_HEIGHT),
                AttributeRule.newIntegerRule(AUXILLIARY_INTERVALS, true),
                new ElementRule(GMRFSkyrideLikelihood.class),
        };

        public XMLSyntaxRule[] getSyntaxRules() {
            return null;
        }

        public String getParserName() {
            return SKYRIDE_FIXED_GRID;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Logger.getLogger("dr.evolmodel").info("Creating GMRF Skyride Fixed Grid Logger");

            GMRFSkyrideLikelihood g = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);

            //Number of intervals
            int intervalNumber = xo.getIntegerAttribute(NUMBER_OF_INTERVALS);
            Logger.getLogger("dr.evomodel").info("Number of Intervals = " + intervalNumber);

            //Grid height
            double gridHeight = xo.getDoubleAttribute(GRID_HEIGHT);
            Logger.getLogger("dr.evomodel").info("Grid Height = " + gridHeight);

            //Filename
            String fileName = null;
            if (xo.hasAttribute(LoggerParser.FILE_NAME)) {
                fileName = xo.getStringAttribute(LoggerParser.FILE_NAME);
            }

            //Log every
            int logEvery = 1;
            if (xo.hasAttribute(LoggerParser.LOG_EVERY)) {
                logEvery = xo.getIntegerAttribute(LoggerParser.LOG_EVERY);
            }

            int auxiliaryIntervals = 1;
            if (xo.hasAttribute(AUXILLIARY_INTERVALS)) {
                auxiliaryIntervals = xo.getIntegerAttribute(AUXILLIARY_INTERVALS);
            }


            PrintWriter pw;

            if (fileName != null) {
                try {
                    File file = new File(fileName);
                    String name = file.getName();
                    String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }
                    pw = new PrintWriter(new FileOutputStream(new File(parent, name)));
                } catch (FileNotFoundException fnfe) {
                    throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
                }
            } else {
                pw = new PrintWriter(System.out);
            }

            LogFormatter lf = new TabDelimitedFormatter(pw);

            MCLogger logger = new GMRFSkyrideFixedGridLogger(lf, logEvery, g, intervalNumber, gridHeight, auxiliaryIntervals);

            //After constructing the logger, setup the columns.
            final BeastVersion version = new BeastVersion();

            logger.setTitle("BEAST " + version.getVersionString() + ", " +
                    version.getBuildString() + "\n" +

                    "Generated " + (new Date()).toString() +
                    "\nFirst value corresponds to coalescent interval closet to sampling time" +
                    "\nLast value corresponds to coalescent interval closet to the root" +
                    "\nGrid Height = " + Double.toString(gridHeight) +
                    "\nNumber of intervals = " + intervalNumber);

            for (int i = 0; i < intervalNumber; i++) {
                logger.addColumn(new LogColumn.Default("V" + (i + 1), null));
            }
            logger.addColumn(new LogColumn.Default("AuxiliaryLikelihood", null));

            return (GMRFSkyrideFixedGridLogger) logger;
        }


    };

}
