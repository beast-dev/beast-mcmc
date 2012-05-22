package dr.evomodel.speciation.ModelAveragingResearch;

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.StringTokenizer;

/**
 * @author Walter Xie
 */
public class MAIndexResearchLog extends LogFileTraces {

    private int[][] indexTranMatrix; // 3 => 6; 4 => 24
    MAIndexResearch ma;

    public MAIndexResearchLog(String name, File file) {
        super(name, file);

        ma = new MAIndexResearch();
        indexTranMatrix = new int[ma.index][ma.index];

        try {
            FileReader reader = new FileReader(file);
            loadTraces(reader);
            reader.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        report();
    }

    private void report() {
        System.out.println("\n-----------------------\n");
        System.out.print(" =>");
        for(int[] p : ma.indexPatterns) {
            System.out.printf("\t%10s", getPatterns(p));
        }
        System.out.print("\n");
        
        int i=0;
        for(int[] row : indexTranMatrix) {
            System.out.print(getPatterns(ma.indexPatterns.get(i)));
            i++;
            for(int col : row) {
                System.out.printf("\t%10d", col);
            }
            System.out.print("\n");
        }
    }

    public void loadTraces(Reader r) throws TraceException, IOException {

        TrimLineReader reader = new LogFileTraces.TrimLineReader(r);

        // Read through to first token
        StringTokenizer tokens = reader.tokenizeLine();

        if (tokens == null) {
            throw new TraceException("Trace file is empty.");
        }

        // read over empty lines
        while (!tokens.hasMoreTokens()) {
            tokens = reader.tokenizeLine();
        }

        // skip the first column which should be the state number
        String token = tokens.nextToken();

        // lines starting with [ are ignored, assuming comments in MrBayes file
        // lines starting with # are ignored, assuming comments in Migrate or BEAST file
        while (token.startsWith("[") || token.startsWith("#")) {

            tokens = reader.tokenizeLine();

            // read over empty lines
            while (!tokens.hasMoreTokens()) {
                tokens = reader.tokenizeLine();
            }

            // read state token and ignore
            token = tokens.nextToken();
        }

        // read label tokens
        String[] labels = new String[tokens.countTokens()];

        for (int i = 0; i < labels.length; i++) {
            labels[i] = tokens.nextToken();
        }

        int[] indexLength = getIndexLength(labels);
        System.out.println("\n-----------------------\n");
        System.out.println("indexLength = " + indexLength[0] + ";   matrixLength = " + indexTranMatrix.length + "\n");

//        int traceCount = getTraceCount();

        boolean firstState = true;

        tokens = reader.tokenizeLine();

        int[] prevValues = new int[0];
        while (tokens != null && tokens.hasMoreTokens()) {

            String stateString = tokens.nextToken();
            int state = 0;

            try {
                try {
                    // Changed this to parseDouble because LAMARC uses scientific notation for the state number
                    state = (int) Double.parseDouble(stateString);
                } catch (NumberFormatException nfe) {
                    throw new TraceException("Unable to parse state number in column 1 (Line " + reader.getLineNumber() + ")");
                }

                if (firstState) {
                    // MrBayes puts 1 as the first state, BEAST puts 0
                    // In order to get the same gap between subsequent samples,
                    // we force this to 0.
                    if (state == 1) state = 0;
                    firstState = false;
                }

            } catch (NumberFormatException nfe) {
                throw new TraceException("State " + state + ":Expected real value in column " + reader.getLineNumber());
            }

            for (int i = 0; i < indexLength[1]; i++) {
                tokens.nextToken();
            }

            int[] values = new int[indexLength[0]];
            for (int i = 0; i < indexLength[0]; i++) {
                try {
                    values[i] = (int) Double.parseDouble(tokens.nextToken());

                } catch (NumberFormatException nfe) {
                    throw new TraceException("State " + state + ": Expected correct number type (Double, Integer or String) in column "
                            + (i + 1) + " (Line " + reader.getLineNumber() + ")");
                }
            }

            if (state > 0) countTranMatrix(prevValues, values);

            prevValues = new int[indexLength[0]];
            System.arraycopy(values, 0, prevValues, 0, values.length);

            tokens = reader.tokenizeLine();
        }
    }

    private void countTranMatrix(int[] prevValues, int[] values) throws TraceException {
        int pre = -1;
        int cur = -1;
        for (int i = 0; i < ma.indexPatterns.size(); i++) {
            int[] p = ma.indexPatterns.get(i);

            if (p.length != prevValues.length) throw new TraceException("MAIndexResearch.indexLen != indexLength in log"); 

            if (ma.isRepeated(p, prevValues)) pre = i;
            if (ma.isRepeated(p, values)) cur = i;
        }

        indexTranMatrix[pre][cur] += 1;  
    }

    private int[] getIndexLength(String[] labels) {
        int[] len = new int[2];        
        for (int i = 0; i < labels.length; i++) {
            if (labels[i].startsWith("integer")) {
                if (len[0] == 0) len[1] = i;
                len[0] += 1;
            }
        }
        return len;
    }

    public String getPatterns(int[] patterns) {
        String p = "";
        for (int pattern : patterns) {
            p += Integer.toString(pattern);
        }
        return p;
    }

    public static void main(String[] args) {
        File file = new File(".\\tanjaCLSTall.log");
        MAIndexResearchLog mal = new MAIndexResearchLog("", file);
    }

}
