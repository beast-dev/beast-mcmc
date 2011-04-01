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
    protected Taxa unionSetTaxonList = new Taxa();
    protected Microsatellite microsatellite;
    protected boolean hasDifferentTaxon = false;

    public MicroSatImporter(BufferedReader reader) {
        this(reader, "\t");
    }

    public MicroSatImporter(BufferedReader reader, String delimiter) {
        this.reader = reader;
        this.delimiter = delimiter;
    }

    public Taxa getUnionSetTaxonList() throws IOException, Importer.ImportException {
        return unionSetTaxonList;
    }

    public Microsatellite getMicrosatellite() {
        return microsatellite;
    }

    public boolean isHasDifferentTaxon() {
        return hasDifferentTaxon;
    }

    public List<Patterns> importPatterns() throws IOException, Importer.ImportException {
        List<Patterns> microsatPatList = new ArrayList<Patterns>();
        List<List<String>> data = new ArrayList<List<String>>(); // 1st List<String> is taxon names
        String[] microsatName = new String[2]; // microsatName[0] is keyword, microsatName[1] is name
        microsatName[1] = "unnamed.microsat";
        String line = reader.readLine();
        while (line.startsWith("#")) { // comments
            if (line.toUpperCase().contains("NAME")) {
                microsatName = line.trim().split("[" + delimiter + " ]+");
                if (microsatName[1] == null || microsatName[1].length() < 1)
                    throw new Importer.ImportException("Improper microsatellite name : " + microsatName[1]);
            }
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

        line = reader.readLine();
        while (line != null) { // read data
            String[] dataLine = line.trim().split("[" + delimiter + " ]+");

            if (dataLine.length != colLen)
                throw new Importer.ImportException("The number of name columns are different with values columns," +
                        "\nplease use only letters or numbers in the name.");
//                + "\ndataLine.length = " + dataLine.length + ", colLen = " + colLen);

            for (int i = 0; i < dataLine.length; i++) {
                data.get(i).add(dataLine[i]);
                if (i > 0) {
                    int v = parseInt(dataLine[i]);
                    if (v != Microsatellite.UNKNOWN_STATE_LENGTH) {
                        if (min > v) min = v;
                        if (max < v) max = v;
                    }
                }
            }

            line = reader.readLine();
        }

        if (max < min) throw new Importer.ImportException("Importing invaild data: max < min !");
//        if (min - 2 < 0) throw new Importer.ImportException("Importing invaild data: min-2 < 0 where min = " + min);
        // The min also = 1 and max should be the longest repeat length + 2.
        microsatellite = new Microsatellite(microsatName[1], 1, max + 2, 1);

        Taxa taxaHaploid = new Taxa();
        for (String name : data.get(0)) {
            Taxon t = new Taxon(name);
            taxaHaploid.addTaxon(t);
        }
//        unionSetTaxonList.addTaxa(taxaHaploid);

        Patterns microsatPat;
        for (int i = 1; i < data.size(); i++) { // create pattern
//            List<Integer> pattern = new ArrayList<Integer>();
            List<Integer> pattern;
            Taxa taxa = new Taxa();

            if ((i + 1 < data.size()) && names[i].equalsIgnoreCase(names[i + 1])) { // diploid: Locus2	Locus2
                Taxa taxaDiploid = new Taxa();
                for (String name : data.get(0)) {
                    Taxon t = new Taxon(names[i] + "_1_" + name);
                    taxaDiploid.addTaxon(t);
                }
                for (String name : data.get(0)) {
                    Taxon t = new Taxon(names[i] + "_2_" + name);
                    taxaDiploid.addTaxon(t);
                }

                if (unionSetTaxonList.containsAny(taxaDiploid))
                    throw new Importer.ImportException("Importing invaild data: duplicate taxon name in this locus : " + names[i]);

                unionSetTaxonList.addTaxa(taxaDiploid);
                hasDifferentTaxon = true;

                pattern = new ArrayList<Integer>();
                String value;
                int size = data.get(i).size();

                for (int v = 0; v < size; v++) {
                    value = data.get(i).get(v);
                    if (!isUnknownChar(value)) {
                        Taxon t = taxaDiploid.getTaxon(v);
                        if (!taxa.contains(t)) {
                            taxa.addTaxon(t);
                            pattern.add(parseInt(value));//microsatellite.getState(value);
                            if (!unionSetTaxonList.contains(t)) {
                                unionSetTaxonList.addTaxon(t);
                                if (i > 1) hasDifferentTaxon = true;
                            }
                        }
                    }
                }
                for (int v = 0; v < data.get(i + 1).size(); v++) {
                    value = data.get(i + 1).get(v);
                    if (!isUnknownChar(value)) {
                        Taxon t = taxaDiploid.getTaxon(v + size);
                        if (!taxa.contains(t)) {
                            taxa.addTaxon(t);
                            pattern.add(parseInt(value));//microsatellite.getState(value);
                            if (!unionSetTaxonList.contains(t)) {
                                unionSetTaxonList.addTaxon(t);
                                if (i > 1) hasDifferentTaxon = true;
                            }
                        }
                    }
                }

                i++;

            } else { // haploid Locus1
                pattern = new ArrayList<Integer>();

                for (int v = 0; v < data.get(i).size(); v++) {
                    String value = data.get(i).get(v);
                    if (!isUnknownChar(value)) {
                        Taxon t = taxaHaploid.getTaxon(v);
                        if (!taxa.contains(t)) {
                            taxa.addTaxon(t);
                            pattern.add(parseInt(value));//microsatellite.getState(value);
                            if (!unionSetTaxonList.contains(t)) {
                                unionSetTaxonList.addTaxon(t);
                                if (i > 1) hasDifferentTaxon = true;
                            }
                        }
                    }
                }

            }
            int[] p = new int[pattern.size()];
            for (int v = 0; v < pattern.size(); v++) {
                p[v] = pattern.get(v);
            }

            microsatPat = new Patterns(microsatellite, taxa);

            microsatPat.addPattern(p);
            microsatPat.setId(names[i]);
            microsatPatList.add(microsatPat);
        }

        return microsatPatList;
    }

    private int parseInt(String s) {
        if (s.charAt(0) == Microsatellite.UNKNOWN_CHARACTER) {
            return Microsatellite.UNKNOWN_STATE_LENGTH; // -1
        } else {
            return Integer.parseInt(s);
        }
    }

    private boolean isUnknownChar(String s) {
        return parseInt(s) == -1; // -1
    }

    /*
id	Locus1	Locus2	Locus2	Locus3	Locus4	Locus4	Locus5	Locus6
T1	5	6	?	20	?	?	?	11
T2	5	6	?	12	?	?	?	12
T3	8	6	4	16	9	9	?	13
T4	12	?	6	1	9	12	?	4
T5	17	?	9	18	7	7	?	5
T6	19	?	5	14	12	12	?	6
     */
}
