/*
 * Citation.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.util;

import java.util.Arrays;

/**
 * @author Alexei Drummond
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class Citation {

    private final Author[] authors;
    private final String title;
    private final int year;
    private final String journal;
    private final String location; // alternative for eJournal
    private final int volume;
    private final int startpage;
    private final int endpage;
    private final Status status;
    private final String DOI;

    public Citation(Author[] authors, Status status) {
        this(authors, null, null, status);
        if (status != Status.IN_PREPARATION) {
            throw new CitationException("Only citations in preparation may not contain titles or journals");
        }
    }

    public Citation(Author[] authors, String title, String journal,
                   Status status) {
        this(authors, title, -1, journal, -1, -1, -1, null, status);
        if (status == Status.PUBLISHED) {
            throw new CitationException("Published citations must have years, volumes and pages");
        }
    }

    public Citation(Author[] authors, String title, int year, String journal, int volume, int startpage, int endpage,
                    Status status) {
        this(authors, title, year, journal, volume, startpage, endpage, null, status);
    }

    public Citation(Author[] authors, String title, int year, String journal, int volume, int startpage, int endpage) {
        this(authors, title, year, journal, volume, startpage, endpage, null, Status.PUBLISHED);
    }

    public Citation(Author[] authors, String title, int year, String journal, int volume, int startpage, int endpage,
                   String DOI) {
        this(authors, title, year, journal, volume, startpage, endpage, DOI, Status.PUBLISHED);
    }

    public Citation(Author[] authors, String title, int year, String journal, int volume, int startpage, int endpage,
                    String DOI, Status status) {
        this.authors = authors;
        this.title = title;
        this.year = year;
        this.journal = journal;
        this.volume = volume;
        this.startpage = startpage;
        this.endpage = endpage;
        this.location = null;
        this.DOI = DOI;
        this.status = status;
    }

    public Citation(Author[] authors, String title, int year, String journal, int volume, String location) {
        this(authors, title, year, journal, volume, location, null);
    }

    public Citation(Author[] authors, String title, int year, String journal, String location) {
        this(authors, title, year, journal, -1, location, null);
    }

    public Citation(Author[] authors, String title, int year, String journal, int volumn, String location,
                    String DOI) {
        this.authors = authors;
        this.title = title;
        this.year = year;
        this.journal = journal;
        this.location = location;
        this.volume = volumn;
        this.startpage = -1;
        this.endpage = -1;
        this.DOI = DOI;
        this.status = Status.PUBLISHED;
    }

    public Citation(Author[] authors, String title, int year, String journal, String location,
                    String DOI) {
        this(authors, title, year, journal, -1, location, DOI);
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append(authors[0].toString());
        for (int i = 1; i < authors.length; i++) {
            builder.append(", ");
            builder.append(authors[i].toString());
        }
        builder.append(" (");
        switch (status) {
            case PUBLISHED: builder.append(year); break;
            default: builder.append(status.getText());
        }
        builder.append(") ");
        if (title != null) {
            builder.append(title);
        }
        if (journal != null) {
            builder.append(". ");
            builder.append(journal);
        }
        if (status == Status.PUBLISHED) {
             builder.append(". ");
            if (location != null) {
                builder.append(location);
            } else {
                builder.append(volume);
                builder.append(", ");
                builder.append(startpage);
                if (endpage > 0) builder.append("-").append(endpage);
            }

            if (DOI != null) {
                builder.append(". DOI:" + DOI);
            }
        }
        return builder.toString();
    }

    public String toHTML() {
        StringBuilder builder = new StringBuilder();
        builder.append("<html>");
        builder.append(authors[0].toString());
        for (int i = 1; i < authors.length; i++) {
            builder.append(", ");
            builder.append(authors[i].toString());
        }
        builder.append(" (").append(year).append(") ");
        builder.append(title).append(". ");
        builder.append("<i>").append(journal).append("</i>");
        if (location != null) {
            builder.append(" ").append(location);
        } else {
            builder.append(" <b>").append(volume).append("</b>:");
            builder.append(startpage);
            if (endpage > 0) builder.append("-").append(endpage);
        }
        if (DOI != null) {
            builder.append(" <a href=\"http://doi.org/").append(DOI).append("\">DOI:").append(DOI).append("</a>");
        }
        builder.append("</html>");

        return builder.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Citation citation = (Citation) o;

        if (year != citation.year) return false;
        if (volume != citation.volume) return false;
        if (startpage != citation.startpage) return false;
        if (endpage != citation.endpage) return false;
        // Probably incorrect - comparing Object[] arrays with Arrays.equals
        if (!Arrays.equals(authors, citation.authors)) return false;
        if (!title.equals(citation.title)) return false;
        if (journal != null ? !journal.equals(citation.journal) : citation.journal != null) return false;
        if (location != null ? !location.equals(citation.location) : citation.location != null) return false;
        if (status != citation.status) return false;
        return DOI != null ? DOI.equals(citation.DOI) : citation.DOI == null;

    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(authors);
        result = 31 * result + title.hashCode();
        result = 31 * result + year;
        result = 31 * result + (journal != null ? journal.hashCode() : 0);
        result = 31 * result + (location != null ? location.hashCode() : 0);
        result = 31 * result + volume;
        result = 31 * result + startpage;
        result = 31 * result + endpage;
        result = 31 * result + status.hashCode();
        result = 31 * result + (DOI != null ? DOI.hashCode() : 0);
        return result;
    }

    public Status getStatus() {
        return status;
    }

    public String getDOI() {
        return DOI;
    }

    public enum Status {
        IN_PREPARATION("in preparation"),
        IN_SUBMISSION("in submission"),
        IN_PRESS("in press"),
        ACCEPTED("accepted"),
        PUBLISHED("published");

        Status(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        private final String text;
    }

    public enum Category {
        FRAMEWORK("Framework"),
        SUBSTITUTION_MODELS("Substitution Models"),
        PRIOR_MODELS("Prior Models"),
        TRAIT_MODELS("Trait Models"),
        DATA_MODELS("Data Models"),
        SPECIES_MODELS("Species Models"),
        COUNTING_PROCESSES("Counting Processes"), // TODO Decide where MarkovJumpsBTL goes (multiple categories?)
        TREE_PRIORS("Tree Density Models"),
        MOLECULAR_CLOCK("Molecular Clock Models"),
        TREE_METRICS("Tree Metrics"),
        MISC("Misc"); // Try to avoid this category


        Category(String text) { this.text = text; }

        public String toString() { return text; }

        private final String text;
    }

    class CitationException extends RuntimeException {
        CitationException(String message) {
            super(message);
        }
    }
}
