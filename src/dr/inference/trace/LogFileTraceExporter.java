package dr.inference.trace;

import dr.util.TabularData;

import java.io.File;
import java.io.IOException;

/**
 * Export trace analysis data such as mean,median,HPD and ESS of trace variables.
 *
 * @author Joseph Heled
 *         Date: 25/10/2007
 */
public class LogFileTraceExporter extends TabularData {
    private final LogFileTraces analysis;
    private final String[] rows = {"mean", "median", "hpdLower", "hpdUpper", "ESS"};
    TraceDistribution[] distributions;

    public LogFileTraceExporter(File file, int burnin) throws TraceException, IOException {

        analysis = new LogFileTraces(file.getCanonicalPath(), file);
        analysis.loadTraces();
        if (burnin >= 0) {
            analysis.setBurnIn(burnin);
        }

        distributions = new TraceDistribution[nColumns()];
    }

    public int nColumns() {
        return analysis.getTraceCount();
    }

    public String columnName(int nColumn) {
        return analysis.getTraceName(nColumn);
    }

    public int nRows() {
        return rows.length;
    }

    public Object data(int nRow, int nColumn) {
        // read on demand
        if (distributions[nColumn] == null) {
            analysis.analyseTrace(nColumn);
            distributions[nColumn] = analysis.getDistributionStatistics(nColumn);
        }

        TraceDistribution distribution = distributions[nColumn];

        switch (nRow) {
            case 0: {
                return distribution.getMean();
            }
            case 1: {
                return distribution.getMedian();
            }
            case 2: {
                return distribution.getLowerHPD();
            }
            case 3: {
                return distribution.getUpperHPD();
            }
            case 4: {
                return distribution.getESS();
            }
        }

        return null;
    }

}
