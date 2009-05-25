package dr.app.beauti.options;

import dr.app.beauti.options.BeautiOptions.*;
import dr.app.beauti.traitspanel.GuessTraitException;
import dr.app.util.Utils;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TraitGuesser {

    public enum GuessType {
        ORDER,
        PREFIX,
        REGEX
    }
    
    public static enum TraitType {
        DISCRETE,
        INTEGER,
        CONTINUOUS
    }
    
    public enum TraitAnalysisType {

    	SPECIES_ANALYSIS ("Species Analysis");

    	TraitAnalysisType(String name) {
    		this.name = name;
    	}

    	public String toString() {
    		return name;
    	}

    	private final String name;
    }
    
    public enum Traits {

    	TRAIT_SPECIES ("species");

    	Traits(String name) {
    		this.name = name;
    	}

    	public String toString() {
    		return name;
    	}

    	private final String name;
    }


    public boolean guessTrait = false;
    public GuessType guessType = GuessType.ORDER;
    public TraitAnalysisType traitAnalysisType = TraitAnalysisType.SPECIES_ANALYSIS;
    public TraitType traitType = TraitType.DISCRETE;
    public boolean fromLast = false;
    public int order = 0;
    public String prefix;
    public String regex;
    public double offset = 0.0;
    public double unlessLessThan = 0.0;
    public double offset2 = 0.0;

    public void guessTrait(BeautiOptions options) {

        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            
            double d = 0.0;

            try {
                switch (guessType) {
                    case ORDER:
                        d = guessTraitFromOrder(options.taxonList.getTaxonId(i), order, fromLast);
                        break;
                    case PREFIX:
                        d = guessTraitFromPrefix(options.taxonList.getTaxonId(i), prefix);
                        break;
                    case REGEX:
                        d = guessTraitFromRegex(options.taxonList.getTaxonId(i), regex);
                        break;
                    default:
                        throw new IllegalArgumentException("unknown GuessType");
                }

            } catch (GuessTraitException gfe) {
                //
            }
            
            options.taxonList.getTaxon(i).setAttribute(Traits.TRAIT_SPECIES.toString(), d);
        }
    }

    public double guessTraitFromOrder(String label, int order, boolean fromLast) throws GuessTraitException {

        String field;

        if (fromLast) {
            int count = 0;
            int i = label.length() - 1;

            char c = label.charAt(i);

            do {
                // first find a part of a number
                while (!Character.isDigit(c) && c != '.') {
                    i--;
                    if (i < 0) break;
                    c = label.charAt(i);
                }

                if (i < 0) throw new GuessTraitException("Missing number field in taxon label, " + label);

                int j = i + 1;

                // now find the beginning of the number
                while (Character.isDigit(c) || c == '.') {
                    i--;
                    if (i < 0) break;
                    c = label.charAt(i);
                }

                field = label.substring(i + 1, j);

                count++;

            } while (count <= order);

        } else {
            int count = 0;
            int i = 0;

            char c = label.charAt(i);

            do {
                // first find a part of a number
                while (!Character.isDigit(c)) {
                    i++;
                    if (i == label.length()) break;
                    c = label.charAt(i);
                }
                int j = i;

                if (i == label.length()) throw new GuessTraitException("Missing number field in taxon label, " + label);

                // now find the beginning of the number
                while (Character.isDigit(c) || c == '.') {
                    i++;
                    if (i == label.length()) break;
                    c = label.charAt(i);
                }

                field = label.substring(j, i);

                count++;

            } while (count <= order);
        }

        return Double.parseDouble(field);
    }

    public double guessTraitFromPrefix(String label, String prefix) throws GuessTraitException {

        int i = label.indexOf(prefix);

        if (i == -1) throw new GuessTraitException("Missing prefix in taxon label, " + label);

        i += prefix.length();
        int j = i;

        // now find the beginning of the number
        char c = label.charAt(i);
        while (i < label.length() - 1 && (Character.isDigit(c) || c == '.')) {
            i++;
            c = label.charAt(i);
        }

        if (i == j) throw new GuessTraitException("Missing field after prefix in taxon label, " + label);

        String field = label.substring(j, i + 1);

        double d;

        try {
            d = Double.parseDouble(field);
        } catch (NumberFormatException nfe) {
            throw new GuessTraitException("Badly formated date in taxon label, " + label);
        }

        return d;
    }

    public double guessTraitFromRegex(String label, String regex) throws GuessTraitException {
        double d;

        try {
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(label);
            if (!matcher.find()) {
                throw new GuessTraitException("Regular expression doesn't find a match in taxon label, " + label);
            }

            if (matcher.groupCount() < 1) {
                throw new GuessTraitException("Date group not defined in regular expression: use parentheses to surround the character group that gives the date");
            }

            d = Double.parseDouble(matcher.group(0));
        } catch (NumberFormatException nfe) {
            throw new GuessTraitException("Badly formated date in taxon label, " + label);
        }

        return d;
    }
}
