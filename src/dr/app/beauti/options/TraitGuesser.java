/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import dr.app.beauti.traitspanel.GuessTraitException;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TraitGuesser {

    private final TraitsOptions traitsOptions;

    public TraitGuesser(String traitName, TraitsOptions.TraitType traitType) {
        this.traitName = traitName;
        this.traitType = traitType;

        if (traitType == TraitsOptions.TraitType.DISCRETE) {
            traitsOptions = new DiscreteTraitOptions(this);
        } else {
            traitsOptions = null; //TODO integer and continuous
        }
    }

    public static enum GuessType {
        SUFFIX,
        PREFIX,
        REGEX
    }

//    private boolean guessTrait = false; // no use ??

    private GuessType guessType = GuessType.SUFFIX;
    private String traitName = TraitsOptions.Traits.TRAIT_SPECIES.toString();
    private TraitsOptions.TraitType traitType = TraitsOptions.TraitType.DISCRETE;

    private int index = 0;
    private String separator;
    private String regex;

    ////////////////////////////////////////////////////////////////

    public TraitsOptions getTraitsOptions() {
        return traitsOptions;
    }
    
    public GuessType getGuessType() {
        return guessType;
    }

    public void setGuessType(GuessType guessType) {
        this.guessType = guessType;
    }

    public String getTraitName() {
        return traitName;
    }

    public void setTraitName(String traitName) {
        this.traitName = traitName;
    }

    public TraitsOptions.TraitType getTraitType() {
        return traitType;
    }

    public void setTraitType(TraitsOptions.TraitType traitType) {
        this.traitType = traitType;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public String getSeparator() {
        return separator;
    }

    public void setSeparator(String separator) {
        this.separator = separator;
    }

    public String getRegex() {
        return regex;
    }

    public void setRegex(String regex) {
        this.regex = regex;
    }

    ////////////////////////////////////////////////////////////////
    public void guessTrait(BeautiOptions options) {

        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {

            String value = null;

            try {
                switch (guessType) {
                    case SUFFIX:
                        value = guessTraitFromSuffix(options.taxonList.getTaxonId(i), separator);
                        break;
                    case PREFIX:
                        value = guessTraitFromPrefix(options.taxonList.getTaxonId(i), separator);
                        break;
                    case REGEX:
                        value = guessTraitFromRegex(options.taxonList.getTaxonId(i), regex);
                        break;
                    default:
                        throw new IllegalArgumentException("unknown GuessType");
                }

            } catch (GuessTraitException gfe) {
                //
            }

            options.taxonList.getTaxon(i).setAttribute(traitName, value);
        }
    }

    public String guessTraitFromSuffix(String label, String seperator) throws GuessTraitException {
        if (seperator.length() < 1 || index < 0) {
            throw new IllegalArgumentException("No seperator or wrong seperator index !");
        }

        int id = -1;
        String t = label;
        for (int i = 0; i <= index; i++) { // i <= index
            id = t.lastIndexOf(seperator);

            if (id < 0) {
                throw new IllegalArgumentException("Can not find seperator in taxon label (" + label + ")\n or invalid seperator (" + seperator + ") !");
            }

            t = t.substring(0, id);
        }

        return label.substring(id + 1); // exclude ith separator
    }

    public String guessTraitFromPrefix(String label, String seperator) throws GuessTraitException {
        if (seperator.length() < 1 || index < 0) {
            throw new IllegalArgumentException("No seperator or wrong seperator index !");
        }

        int id;
        int idSum = 0;
        String t = label;
        for (int i = 0; i <= index; i++) { // i <= index
            id = t.indexOf(seperator);

            if (id < 0) {
            	if (i == 0) {
            		throw new IllegalArgumentException("Can not find seperator in taxon label (" + label + ")\n or invalid seperator (" + seperator + ") !");
            	} else {
            		return label + seperator + traitName;
            	}
            }

            t = t.substring(id + 1);
            if (i == 0) {
                idSum = idSum + id;
            } else {
                idSum = idSum + id + 1;
            }
        }

        return label.substring(0, idSum); // exclude ith separator
    }

    public String guessTraitFromRegex(String label, String regex) throws GuessTraitException {
        String t;

        try {
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(label);
            if (!matcher.find()) {
                throw new GuessTraitException("Regular expression doesn't find a match in taxon label, " + label);
            }

            if (matcher.groupCount() < 1) {
                throw new GuessTraitException("Trait value group not defined in regular expression");
            }

            t = matcher.group(0); // TODO: not working?
        } catch (NumberFormatException nfe) {
            throw new GuessTraitException("Badly formated trait value in taxon label, " + label);
        }

        return t;
    }
}
