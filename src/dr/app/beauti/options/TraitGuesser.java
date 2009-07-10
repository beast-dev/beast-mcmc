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

    public enum GuessType {
        SUFFIX,
        PREFIX,
        REGEX
    }

    public static enum TraitType {
        DISCRETE,
        INTEGER,
        CONTINUOUS
    }

    public enum Traits {

        TRAIT_SPECIES("species");

        Traits(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }

    public boolean guessTrait = false; // no use ??
    public GuessType guessType = GuessType.SUFFIX;
    public Traits traitAnalysisType = Traits.TRAIT_SPECIES;
    public TraitType traitType = TraitType.DISCRETE;

    public int index = 0;
    public String separator;
    public String regex;

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

            options.taxonList.getTaxon(i).setAttribute(Traits.TRAIT_SPECIES.toString(), value);
        }
    }

    public String guessTraitFromSuffix(String label, String seperator) throws GuessTraitException {
        if (seperator.length() < 1 || index < 0) {
            throw new IllegalArgumentException("Invalid seperator");
        }

        int id = -1;
        String t = label;
        for (int i = 0; i <= index; i++) { // i <= index
            id = t.lastIndexOf(seperator);

            if (id < 0) {
                throw new IllegalArgumentException("Can not find seperator in taxon label or invalid seperator index");
            }

            t = t.substring(0, id);
        }

        return label.substring(id + 1); // exclude ith separator
    }

    public String guessTraitFromPrefix(String label, String seperator) throws GuessTraitException {
        if (seperator.length() < 1 || index < 0) {
            throw new IllegalArgumentException("Invalid seperator");
        }

        int id;
        int idSum = 0;
        String t = label;
        for (int i = 0; i <= index; i++) { // i <= index
            id = t.indexOf(seperator);

            if (id < 0) {
                throw new IllegalArgumentException("Can not find seperator in taxon label or invalid seperator index");
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
