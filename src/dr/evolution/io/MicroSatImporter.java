package dr.evolution.io;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class MicroSatImporter implements PatternImporter {
    protected final BufferedReader reader;
    protected String delimiter;

    public MicroSatImporter(BufferedReader reader) {
        this(reader, "\t");
    }

    public MicroSatImporter(BufferedReader reader, String delimiter) {
        this.reader = reader;
        this.delimiter = delimiter;
    }

    public List<Patterns> importPatterns() throws IOException, Importer.ImportException {
        List<Patterns> pList = new ArrayList<Patterns>();
        List<List<String>> data = new ArrayList<List<String>>(); // 1st List<String> is taxon names

        String line = reader.readLine();
        while (line.startsWith("#") || line.startsWith("[")) { // comments
            line = reader.readLine();
        }
        // read locus (microsat pattern) names in the 1st row after comments, where 1st element is id
        String[] names = line.trim().split("[" + delimiter + " ]+"); // trim trailing whitespace ?
        int colLen = names.length; // for validation

        if (colLen < 2) throw new Importer.ImportException("Import file must have more than 1 columns : " + colLen);

        for (int i = 0; i < colLen; i++) { // init data
            List<String> l = new ArrayList<String>();
            data.add(l);
        }

        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;

        while (line != null) { // read data
            String[] dataLine = line.trim().split("[" + delimiter + " ]+");

            if (dataLine.length != colLen)
                throw new Importer.ImportException("It contains different number of values in Taxon " + dataLine[0]);

            for (int i = 0; i < dataLine.length; i++) {
                data.get(i).add(dataLine[i]);
                int v = parseInt(dataLine[i]);
                if (i > 0 && v != Microsatellite.UNKNOWN_STATE_LENGTH) {
                    if (min > v) min = v;
                    if (max < v) max = v;
                }
            }
        }

        if (max < min) throw new Importer.ImportException("The max < min !");
        if (min - 2 < 0) throw new Importer.ImportException("min-2 < 0 where min = " + min);
        // The min should be the shortest repeat length - 2 and max should be the longest repeat length - 2.
        Microsatellite microsatellite = new Microsatellite(min - 2, max - 2, 1);

        for (int i = 1; i < data.size(); i++) { // create pattern
//            List<Integer> pattern = new ArrayList<Integer>();
            int[] pattern;
            Taxa taxa = new Taxa();

            if ((i + 1 < data.size()) && names[i].equalsIgnoreCase(names[i + 1])) { // e.g. Locus2	Locus2
                pattern = new int[data.get(i).size() + data.get(i + 1).size()]; // todo Jessie ?

                for (int v = 0; v < data.get(i).size(); v++) {
                    String value = data.get(i).get(v);
//            if (parseInt(value) >= 0) { // todo getState handling unused taxon?
                    Taxon t = new Taxon(data.get(0).get(v));
                    taxa.addTaxon(t);
                    pattern[v] = microsatellite.getState(value);
//            }
                }
                for (int v = 0; v < data.get(i + 1).size(); v++) {
                    String value = data.get(i + 1).get(v);
//            if (parseInt(value) >= 0) { // todo getState handling unused taxon?
                    pattern[v] = microsatellite.getState(value);
//            }
                }

                i++;

            } else {
                pattern = new int[data.get(i).size()];

                for (int v = 0; v < data.get(i).size(); v++) {
                    String value = data.get(i).get(v);
//            if (parseInt(value) >= 0) { // todo getState handling unused taxon?
                    Taxon t = new Taxon(data.get(0).get(v));
                    taxa.addTaxon(t);
                    pattern[v] = microsatellite.getState(value);
//            }
                }

            }

            Patterns microsatPat = new Patterns(microsatellite, taxa);
            microsatPat.addPattern(pattern);
            microsatPat.setId(names[i]);
            pList.add(microsatPat);
        }

        return pList;
    }

    private int parseInt(String s) {
        if (s.charAt(0) == Microsatellite.UNKNOWN_CHARACTER) {
            return Microsatellite.UNKNOWN_STATE_LENGTH; // -1
        } else {
            return Integer.parseInt(s);
        }
    }

    /*
id	Locus1	Locus2	Locus2	Locus3	Locus4	Locus4	Locus5	Locus6
T1	5	6	?	20	?	?	?	1
T2	5	6	?	2	?	?	?	2
T3	8	6	4	16	9	9	?	3
T4	12	?	6	1	1	1	?	4
T5	17	?	9	18	7	7	?	5
T6	19	?	5	14	2	2	?	6
     */
}
