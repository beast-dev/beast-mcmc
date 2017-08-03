package dr.evolution.sequence;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Marc A. Suchard
 */
public class UncertainSequence extends Sequence {

    public static char openDelimitor = '{';
    public static char closeDelimitor = '}';
    public static String separateDelimitor = "|";
    public static String traitDelimitor = ":";
    public static char defaultAmbiguousChar = '?';
    
    public char getChar(int index) {

        checkParsed();
        UncertainCharacterList list = characters.get(index);
        if (list.size() == 1) {
            return list.get(0).getCharacter();
        } else {
            return defaultAmbiguousChar;
        }
    }

    public int getLength() {
        checkParsed();
        return characters.size();
    }

    @Override
    public int getState(int index) {

        checkParsed();
        UncertainCharacterList list = characters.get(index);
        if (list.size() == 1) {
            return dataType.getState(list.get(0).getCharacter());
        } else {
            return dataType.getState(defaultAmbiguousChar);
        }
    }


    @Override
    public void setState(int index, int state) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getInvalidChar() {

        checkParsed();

        final char[] validChars = dataType.getValidChars();
        if (validChars != null) {
            String validSet = new String(validChars);

            int index = 0;
            for (UncertainCharacterList list : characters) {
                for (UncertainCharacter testChar : list) {
                    if (!testChar.isValidCharacter(validSet)) {
                        return index;
                    }
                }
                ++index;
            }
        }
        return -1;
    }

    public static boolean containsWeights(String sequenceString) {
        return (sequenceString.indexOf(openDelimitor) != -1) &&
                (sequenceString.indexOf(closeDelimitor) != -1);
    }

    private void checkParsed() {
        if (!isParsed) {
            parseSequenceString();
        }
    }

    private void parseSequenceString() {  // TODO Should throw exceptions if unable to parse information

        characters = new ArrayList<UncertainCharacterList>();

        final String string = sequenceString.toString();
        int bufferIndex = 0;

        while (bufferIndex < string.length()) {
            char current = string.charAt(bufferIndex);

            int end = -1;
            if (current == openDelimitor) {
                end = string.indexOf(closeDelimitor, bufferIndex);
            }

            UncertainCharacterList list = new UncertainCharacterList();
            if (end != -1) {

                String subString = string.substring(bufferIndex + 1, end);
                StringTokenizer weights = new StringTokenizer(subString, separateDelimitor);
                while (weights.hasMoreTokens()) {
                    String element = weights.nextToken();

                    StringTokenizer traits = new StringTokenizer(element, traitDelimitor);
                    String charString = traits.nextToken();
                    String weightString = traits.nextToken();

                    list.add(new UncertainCharacter(charString.charAt(0), Double.valueOf(weightString)));
                }
                bufferIndex = end + 1;

            } else {

                list.add(new UncertainCharacter(current));
                ++bufferIndex;

            }
            characters.add(list);
        }
        isParsed = true;
    }

    public double[] getUncertainPattern(int siteIndex) {
        double[] pattern = new double[dataType.getStateCount()];
        UncertainCharacterList list = characters.get(siteIndex);
        for (UncertainCharacter uc : list) {
            int[] states = dataType.getStates(dataType.getState(uc.getCharacter()));
            for (int state : states) {
                pattern[state] = uc.getWeight();
            }
        }

        return pattern;
    }

    private class UncertainCharacter {
        char character;
        double weight;

        UncertainCharacter(char character) {
            this.character = character;
            this.weight = 1.0;
        }

        UncertainCharacter(char character, double weight) {
            this.character = character;
            this.weight = weight;
        }

        public boolean isValidCharacter(String validSet) {
            return (validSet.indexOf(character) >= 0);
        }

        public char getCharacter() {
            return character;
        }

        public double getWeight() {
            return weight;
        }
    }

    private class UncertainCharacterList extends ArrayList<UncertainCharacter> { }

    private List<UncertainCharacterList> characters;

    private boolean isParsed = false;
}
