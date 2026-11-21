/*
 * FixedEpochProviderParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
 * http://beast.community/about
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodelxml.epochs;

import dr.evolution.util.Date;
import dr.evolution.util.TaxonList;
import dr.evolution.util.TimeScale;
import dr.evolution.util.Units;
import dr.evomodel.epochs.FixedEpochProvider;
import dr.xml.*;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.TimeZone;
import java.util.logging.Logger;

/**
 * Epochs with fixed transition times, given either as a number of epochs of the same duration:
 * <pre>
 *     &lt;fixedEpochs id="epochs" epochCount="12" duration="1" units="months"/&gt;
 * </pre>
 * or as the transition times themselves:
 * <pre>
 *     &lt;fixedEpochs id="epochs" transitionTimes="1 2 4 8 26" units="weeks"/&gt;
 * </pre>
 * The durations or times are in "units" (days, weeks, months or years; years if it is not given) and are
 * converted to years, the units of time of the tree.
 *
 * The transition times are normally counted back from the most recent tip. With the taxa (or a tree) that
 * have dates inside the element, "origin" gives a date, on the scale of the taxon dates, from which they
 * are counted back instead:
 * <pre>
 *     &lt;fixedEpochs id="epochs" epochCount="12" duration="1" units="months" origin="2020.0"&gt;
 *         &lt;taxa idref="taxa"/&gt;
 *     &lt;/fixedEpochs&gt;
 * </pre>
 * puts the transitions one month, two months ... before the start of 2020, whatever the date of the most
 * recent tip. The origin can also be a calendar date, "2020-01-01", when the dates of the taxa are in years.
 *
 * The transitions can be given as calendar dates too, in any order, again with the taxa or a tree whose
 * dates are in years:
 * <pre>
 *     &lt;fixedEpochs id="epochs" transitionTimes="2020-03-23 2020-07-04 2020-11-05"&gt;
 *         &lt;taxa idref="taxa"/&gt;
 *     &lt;/fixedEpochs&gt;
 * </pre>
 * A calendar date is the start of that day and is converted to a decimal year as BEAUti converts the dates
 * of tips.
 *
 * @author Andrew Rambaut
 */
public class FixedEpochProviderParser extends AbstractXMLObjectParser {

    public static final String FIXED_EPOCHS = FixedEpochProvider.FIXED_EPOCHS;
    public static final String EPOCH_COUNT = "epochCount";
    public static final String DURATION = "duration";
    public static final String TRANSITION_TIMES = FixedEpochProvider.TRANSITION_TIMES;
    public static final String UNITS = "units";
    public static final String ORIGIN = "origin";

    private static final String[] TIME_UNIT_NAMES = {"days", "weeks", "months", "years"};
    private static final Units.Type[] TIME_UNIT_TYPES =
            {Units.Type.DAYS, Units.Type.WEEKS, Units.Type.MONTHS, Units.Type.YEARS};

    // the time of a tree is in years
    private static final Units.Type TREE_UNITS = Units.Type.YEARS;

    public String getParserName() {
        return FIXED_EPOCHS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final boolean hasEqualDurations = xo.hasAttribute(EPOCH_COUNT) || xo.hasAttribute(DURATION);
        final boolean hasTransitionTimes = xo.hasAttribute(TRANSITION_TIMES);
        if (hasEqualDurations == hasTransitionTimes) {
            throw new XMLParseException("The " + FIXED_EPOCHS + " element" + idText(xo) + " needs either " +
                    EPOCH_COUNT + " and " + DURATION + ", or " + TRANSITION_TIMES + " (and not both)");
        }

        // the transition times are either all numbers or all calendar dates
        final String[] calendarDates = hasTransitionTimes ? getCalendarDates(xo) : null;
        final boolean hasCalendarDates = calendarDates != null;

        double[] transitionTimes;
        try {
            if (hasEqualDurations) {
                if (!xo.hasAttribute(EPOCH_COUNT) || !xo.hasAttribute(DURATION)) {
                    throw new XMLParseException("The " + FIXED_EPOCHS + " element" + idText(xo) + " needs both " +
                            EPOCH_COUNT + " and " + DURATION);
                }
                transitionTimes = FixedEpochProvider.equalDurationTransitionTimes(
                        xo.getIntegerAttribute(EPOCH_COUNT), xo.getDoubleAttribute(DURATION));
            } else if (!hasCalendarDates) {
                transitionTimes = xo.getDoubleArrayAttribute(TRANSITION_TIMES);
            } else {
                transitionTimes = null;
            }
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ": " + iae.getMessage());
        }

        final Logger logger = Logger.getLogger("dr.evomodel");
        final StringBuilder notes = new StringBuilder();

        if (hasCalendarDates) {
            if (xo.hasAttribute(UNITS) || xo.hasAttribute(ORIGIN)) {
                throw new XMLParseException("The " + FIXED_EPOCHS + " element" + idText(xo) + " gives its " +
                        TRANSITION_TIMES + " as calendar dates, which need no " + UNITS + " or " + ORIGIN);
            }
            // a date is the time before the most recent tip that it falls
            final Date mostRecent = getMostRecentCalendarDate(xo, "calendar dates for " + TRANSITION_TIMES);
            transitionTimes = new double[calendarDates.length];
            for (int i = 0; i < calendarDates.length; i++) {
                transitionTimes[i] = mostRecent.getTimeValue() - parseCalendarDate(calendarDates[i], xo);
            }
            Arrays.sort(transitionTimes);
            if (transitionTimes.length > 0 && !(transitionTimes[0] > 0.0)) {
                throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", the latest of the " +
                        TRANSITION_TIMES + " is not before the most recent tip (" + mostRecent.getTimeValue() + ")");
            }
            String[] sortedDates = calendarDates.clone();
            Arrays.sort(sortedDates, java.util.Collections.reverseOrder());
            notes.append("  transitions at ").append(String.join(", ", sortedDates))
                    .append("; the most recent tip is at ").append(mostRecent.getTimeValue()).append("\n");
        } else {
            if (xo.hasAttribute(UNITS)) {
                transitionTimes = convert(transitionTimes, parseTimeUnits(xo.getStringAttribute(UNITS), xo), TREE_UNITS);
            }

            if (xo.hasAttribute(ORIGIN)) {
                // the transition times are counted back from the origin, not from the most recent tip
                final String originAttribute = xo.getStringAttribute(ORIGIN).trim();
                final Date mostRecent;
                final double origin;
                if (isCalendarDate(originAttribute)) {
                    mostRecent = getMostRecentCalendarDate(xo, "a calendar date for the " + ORIGIN);
                    origin = parseCalendarDate(originAttribute, xo);
                } else {
                    mostRecent = getMostRecentDate(xo);
                    try {
                        origin = Double.parseDouble(originAttribute);
                    } catch (NumberFormatException nfe) {
                        throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", " + ORIGIN +
                                "=\"" + originAttribute + "\" is neither a number nor a calendar date (yyyy-mm-dd)");
                    }
                }
                // dates run forwards in time (2020.5) or backwards (an age of 3.5)
                final double timeBeforeMostRecent = mostRecent.isBackwards() ?
                        origin - mostRecent.getTimeValue() : mostRecent.getTimeValue() - origin;
                final double originHeight =
                        TimeScale.convertTimeUnits(timeBeforeMostRecent, mostRecent.getUnits(), TREE_UNITS);
                for (int i = 0; i < transitionTimes.length; i++) {
                    transitionTimes[i] += originHeight;
                }
                if (transitionTimes.length > 0 && !(transitionTimes[0] > 0.0)) {
                    throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", the first " +
                            "transition between epochs is not before the most recent tip (" +
                            mostRecent.getTimeValue() + "): with an " + ORIGIN + " of " + originAttribute +
                            " it is " + (-transitionTimes[0]) + " years after it.");
                }
                notes.append("  counted back from ").append(originAttribute)
                        .append(isCalendarDate(originAttribute) ? " (" + origin + ")" : "").append(", ")
                        .append(Math.abs(originHeight)).append(" years ")
                        .append(originHeight >= 0.0 ? "before" : "after")
                        .append(" the most recent tip (").append(mostRecent.getTimeValue()).append(")\n");
            }
        }

        final FixedEpochProvider epochs;
        try {
            epochs = new FixedEpochProvider(transitionTimes);
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ": " + iae.getMessage());
        }

        logger.info("\nFixed epochs" + idText(xo) + ": " + epochs.getEpochCount() +
                (epochs.getEpochCount() == 1 ? " epoch" : " epochs"));
        if (transitionTimes.length > 0) {
            logger.info(notes + "  transition times (years before the most recent tip): " +
                    Arrays.toString(transitionTimes));
        }

        return epochs;
    }

    /**
     * @return the times converted from one unit of time to another
     */
    public static double[] convert(double[] times, Units.Type fromUnits, Units.Type toUnits) {
        double[] converted = new double[times.length];
        for (int i = 0; i < times.length; i++) {
            converted[i] = TimeScale.convertTimeUnits(times[i], fromUnits, toUnits);
        }
        return converted;
    }

    /**
     * @return the tokens of the transitionTimes attribute if they are calendar dates, or null if they are
     * numbers
     */
    private static String[] getCalendarDates(XMLObject xo) throws XMLParseException {
        final String value = xo.getStringAttribute(TRANSITION_TIMES).trim();
        if (value.isEmpty()) {
            return null;
        }
        final String[] tokens = value.split("[\\s,]+");
        int dateCount = 0;
        for (String token : tokens) {
            if (isCalendarDate(token)) {
                dateCount++;
            }
        }
        if (dateCount == 0) {
            return null;
        }
        if (dateCount < tokens.length) {
            throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", the " +
                    TRANSITION_TIMES + " must be either all numbers or all calendar dates (yyyy-mm-dd)");
        }
        return tokens;
    }

    private static boolean isCalendarDate(String value) {
        return value.matches("\\d\\d\\d\\d-\\d\\d-\\d\\d");
    }

    /**
     * @return the calendar date as a decimal year, by the conversion BEAUti uses for the dates of tips
     */
    private static double parseCalendarDate(String value, XMLObject xo) throws XMLParseException {
        final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
        format.setLenient(false);
        format.setTimeZone(TimeZone.getTimeZone("GMT"));
        try {
            return new Date(format.parse(value)).getTimeValue();
        } catch (ParseException pe) {
            throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", \"" + value +
                    "\" is not a calendar date");
        }
    }

    /**
     * @return the date of the most recent of the taxa, which for calendar dates to be placed among them
     * must be in years running forwards
     */
    private static Date getMostRecentCalendarDate(XMLObject xo, String purpose) throws XMLParseException {
        if (xo.getChild(TaxonList.class) == null) {
            throw new XMLParseException("The " + FIXED_EPOCHS + " element" + idText(xo) + " has " + purpose +
                    " but not the taxa or tree whose dates place it. Put them inside the element.");
        }
        final Date mostRecent = getMostRecentDate(xo);
        if (mostRecent.getUnits() != Units.Type.YEARS || mostRecent.isBackwards()) {
            throw new XMLParseException("The " + FIXED_EPOCHS + " element" + idText(xo) + " has " + purpose +
                    ", which needs the dates of the taxa to be in years running forwards (such as 2020.5)");
        }
        return mostRecent;
    }

    /**
     * @return the date of the most recent of the taxa or tree inside the element. The dates must all be in
     * the same units and run in the same direction, so that an origin can be on their scale.
     */
    private static Date getMostRecentDate(XMLObject xo) throws XMLParseException {
        final TaxonList taxa = (TaxonList) xo.getChild(TaxonList.class);
        if (taxa == null) {
            throw new XMLParseException("The " + FIXED_EPOCHS + " element" + idText(xo) + " has an " + ORIGIN +
                    " but not the taxa or tree whose dates it is on the scale of. Put them inside the element.");
        }
        Date mostRecent = null;
        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            final Date date = taxa.getTaxon(i).getDate();
            if (date == null) {
                continue;
            }
            if (mostRecent == null) {
                mostRecent = date;
            } else {
                if (date.getUnits() != mostRecent.getUnits() || date.isBackwards() != mostRecent.isBackwards()) {
                    throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", the dates of " +
                            "the taxa are not all in the same units and direction, so dates cannot be placed " +
                            "among them");
                }
                // going backwards the most recent has the smallest value, going forwards the largest
                if (date.isBackwards() ? date.getTimeValue() < mostRecent.getTimeValue() :
                        date.getTimeValue() > mostRecent.getTimeValue()) {
                    mostRecent = date;
                }
            }
        }
        if (mostRecent == null) {
            throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", the taxa have no " +
                    "dates, so dates cannot be placed among them");
        }
        if (!Arrays.asList(TIME_UNIT_TYPES).contains(mostRecent.getUnits())) {
            throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", the dates of the " +
                    "taxa are not in units of time, so dates cannot be placed among them");
        }
        return mostRecent;
    }

    private static Units.Type parseTimeUnits(String name, XMLObject xo) throws XMLParseException {
        for (int i = 0; i < TIME_UNIT_NAMES.length; i++) {
            if (TIME_UNIT_NAMES[i].equalsIgnoreCase(name.trim())) {
                return TIME_UNIT_TYPES[i];
            }
        }
        throw new XMLParseException("In the " + FIXED_EPOCHS + " element" + idText(xo) + ", " + UNITS + "=\"" +
                name + "\" is not a unit of time (days, weeks, months or years)");
    }

    private static String idText(XMLObject xo) {
        try {
            return xo.hasId() ? " '" + xo.getId() + "'" : "";
        } catch (XMLParseException e) {
            return "";
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Epochs of time with fixed transition times, given as a number of epochs of the same duration " +
                "or as the transition times.";
    }

    public Class getReturnType() {
        return FixedEpochProvider.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(EPOCH_COUNT, true,
                    "The number of epochs, counting the last, which has no end"),
            AttributeRule.newDoubleRule(DURATION, true,
                    "The duration of each epoch"),
            new StringAttributeRule(TRANSITION_TIMES,
                    "The times at which one epoch gives way to the next, in increasing order, or their " +
                            "calendar dates (yyyy-mm-dd) in any order", true),
            new StringAttributeRule(UNITS,
                    "The units of the duration or transition times (years if not given)",
                    TIME_UNIT_NAMES, true),
            new StringAttributeRule(ORIGIN,
                    "A date, on the scale of the dates of the taxa or as a calendar date (yyyy-mm-dd), that the " +
                            "transition times are counted back from (by default, the date of the most recent tip)",
                    true),
            new ElementRule(TaxonList.class,
                    "The taxa, or a tree, with the dates that an origin or calendar dates are placed among", 0, 1)
    };
}
