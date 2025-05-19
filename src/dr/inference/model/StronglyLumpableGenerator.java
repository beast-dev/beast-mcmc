package dr.inference.model;

import dr.evolution.datatype.DataType;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collections;

public class StronglyLumpableGenerator {
    private final String xmlOutput;
    private final java.util.List<String> parameterRateIds = new java.util.ArrayList<>();
    private final java.util.List<String> parameterProportionIds = new java.util.ArrayList<>();


    public StronglyLumpableGenerator(StateSet[] stateSets, DataType dataType) {
        this.xmlOutput = generateXML(stateSets, dataType);
        saveToFile("generated_strongly_lumpable_block.xml");
    }
    private String generateXML(StateSet[] stateSets, DataType dataType) {
        StringBuilder xml = new StringBuilder();
        xml.append("<stronglyLumpableCtmcRates id=\"lumpableRates\" normalize=\"true\">\n");

//        // Collect all unique state names
//        Set<String> uniqueStates = new HashSet<>();
//        for (StateSet stateSet : stateSets) {
//            for (int stateId : stateSet.states()) {
//                String stateName = dataType.getCode(stateId);
//                uniqueStates.add(stateName);
//            }
//        }

        // Add the generalDataType block
        String dataTypeName = dataType.getDescription();
        xml.append("\t<generalDataType idref=\"").append(dataTypeName).append("\"/>\n");
//        for (String state : uniqueStates) {
//            xml.append("\t\t<state code=\"").append(state).append("\"/>\n");
//        }
//        xml.append("\t</generalDataType>\n");



        // across lump rate
        xml.append("\t<rates>\n");
        int acrossRateValueSize = stateSets.length * (stateSets.length - 1);
        String acrossRateValues = new String(new char[acrossRateValueSize]).replace("\0", "1 ").trim();
        xml.append("\t\t<parameter id=\"lump.across.rates\"  value=\"").append(acrossRateValues).append("\" lower=\"0\"/>\n");
        parameterRateIds.add("lump.across.rates");
        xml.append("\t</rates>\n");


//        for (StateSet stateSet : stateSets) {
//            // Modify stateSet block inside each lump
//            xml.append("\t<stateSet id=\"").append(stateSet.getId()).append("\">\n");
//            xml.append("\t\t<generalDataType idref=\"loc.dataType\"/>\n");
//            for (String state : stateSet.statesNameArray()) {
//                xml.append("\t\t<state code=\"").append(state).append("\"/>\n");
//            }
//            xml.append("\t</stateSet>\n");
//        }

        //lump
        for (int i = 0; i < stateSets.length; i++) {
            String lumpId = "L" + (i + 1);
            int numStates = stateSets[i].states().length;
            int rateValueSize = numStates * (numStates - 1);
            String rateValues = new String(new char[rateValueSize]).replace("\0", "1 ").trim();

            xml.append("\t<lump id=\"").append(lumpId).append("\">\n");
            xml.append("\t\t<stateSet idref=\"").append(stateSets[i].getId()).append("\"/>\n");

            // Handle <rates> block
            if (rateValues.isEmpty()) {
                xml.append("\t\t<!-- \n");
                xml.append("\t\t<rates>\n");
                xml.append("\t\t\t<parameter id=\"").append(lumpId).append(".within.rates\" value=\"").append(rateValues).append("\" lower=\"0\"/>\n");
                xml.append("\t\t</rates>\n");
                xml.append("\t\t -->\n");
            } else {
                xml.append("\t\t<rates>\n");
                //xml.append("\t\t\t<parameter id=\"").append(lumpId).append(".within.rates\" value=\"").append(rateValues).append("\"/>\n");
                String withinRateId = lumpId + ".within.rates";
                parameterRateIds.add(withinRateId);
                xml.append("\t\t\t<parameter id=\"").append(withinRateId).append("\" value=\"").append(rateValues).append("\" lower=\"0\"/>\n");

                xml.append("\t\t</rates>\n");
            }

            // Handle <proportion> blocks

            for (int stateId : stateSets[i].states()) { // Retrieve state IDs from stateSet
                String stateName = dataType.getCode(stateId);

                for (int j = 0; j < stateSets.length; j++) {
                    if (j != i) {
                        int proportionValueSize = stateSets[j].states().length;
                        // proportions should sum up to 1, so just set 1/proportionSize as default
                        double proportion = 1.0 / proportionValueSize;
                        DecimalFormat df = new DecimalFormat("0.00");
                        String proportionStr = df.format(proportion);
                        String proportionValues = String.join(" ", Collections.nCopies(proportionValueSize, proportionStr));

                        if (proportionValues.isEmpty()) {
                            xml.append("\t\t<!-- \n");
                            xml.append("\t\t<proportions>\n");
                            xml.append("\t\t\t<state code=\"").append(stateName).append("\"/>\n");
                            xml.append("\t\t\t<stateSet idref=\"").append(stateSets[j].getId()).append("\"/>\n");
                            xml.append("\t\t\t<parameter id=\"").append(stateName).append(".").append(stateSets[j].getId()).append(".proportions\" value=\"").append(proportionValues).append("\" lower=\"0\" upper=\"1\"/>\n");
                            xml.append("\t\t</proportions>\n");
                            xml.append("\t\t -->\n");
                        } else {
                            xml.append("\t\t<proportions>\n");
                            xml.append("\t\t\t<state code=\"").append(stateName).append("\"/>\n");
                            xml.append("\t\t\t<stateSet idref=\"").append(stateSets[j].getId()).append("\"/>\n");
                            //xml.append("\t\t\t<parameter id=\"").append(stateName).append(".").append(stateSets[j].getId()).append(".proportions\" value=\"").append(proportionValues).append("\"/>\n");
                            String proportionId = stateName + "." + stateSets[j].getId() + ".proportions";
                            parameterProportionIds.add(proportionId);
                            xml.append("\t\t\t<parameter id=\"").append(proportionId).append("\" value=\"").append(proportionValues).append("\" lower=\"0\" upper=\"1\"/>\n");
                            xml.append("\t\t</proportions>\n");
                        }
                    }
                }
            }

            xml.append("\t</lump>\n");
        }

        xml.append("</stronglyLumpableCtmcRates>");

        // operator block
        xml.append("\n<operators>\n");
        for (String id : parameterRateIds) {
            xml.append("\t<scaleOperator scaleFactor=\"0.75\" weight=\"1\">\n");
            xml.append("\t\t<parameter idref=\"").append(id).append("\"/>\n");
            xml.append("\t</scaleOperator>\n");
        }
        for (String id : parameterProportionIds) {
            xml.append("\t<deltaExchange delta=\"0.75\" weight=\"1\">\n");
            xml.append("\t\t<parameter idref=\"").append(id).append("\"/>\n");
            xml.append("\t</deltaExchange>\n");
        }
        xml.append("</operators>\n");

        // Add prior block
        xml.append("\n<prior>\n");
        for (String id : parameterRateIds) {
            xml.append("\t<gammaPrior shape=\"0.05\" scale=\"10.0\" offset=\"0.0\">\n");
            xml.append("\t\t<parameter idref=\"").append(id).append("\"/>\n");
            xml.append("\t</gammaPrior>\n");
        }
        for (String id : parameterProportionIds) {
            xml.append("\t<uniformPrior lower=\"0\" upper=\"1\">\n");
            xml.append("\t\t<parameter idref=\"").append(id).append("\"/>\n");
            xml.append("\t</uniformPrior>\n");
        }
        xml.append("</prior>\n");

        // Add log block
        xml.append("\n<log id=\"screenLog\">\n");
        for (String id : parameterRateIds) {
            xml.append("\t<column label=\"").append(id).append("\" dp=\"4\" width=\"12\">\n");
            xml.append("\t\t<likelihood idref=\"").append(id).append("\"/>\n");
            xml.append("\t</column>\n");
        }
        for (String id : parameterProportionIds) {
            xml.append("\t<column label=\"").append(id).append("\" dp=\"4\" width=\"12\">\n");
            xml.append("\t\t<likelihood idref=\"").append(id).append("\"/>\n");
            xml.append("\t</column>\n");
        }
        xml.append("</log>\n");

        // Add log block for fileLog
        xml.append("\n<log id=\"fileLog\">\n");
        for (String id : parameterRateIds) {
            xml.append("\t<parameter idref=\"").append(id).append("\"/>\n");
        }
        for (String id : parameterProportionIds) {
            xml.append("\t<parameter idref=\"").append(id).append("\"/>\n");
        }
        xml.append("</log>\n");

        return xml.toString();
    }

    public String getXmlOutput() {
        return xmlOutput;
    }


    private void saveToFile(String filename) {
        try (FileWriter fileWriter = new FileWriter(filename)) {
            fileWriter.write(xmlOutput);
            System.out.println("XML file saved as: " + filename);
        } catch (IOException e) {
            System.err.println("Error writing XML to file: " + e.getMessage());
        }
    }
}

