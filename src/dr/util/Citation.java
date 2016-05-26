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

/**
 * @author Alexei Drummond
 * @author Marc A. Suchard
 */
public class Citation {

    Author[] authors;
    String title;
    int year;
    String journal;
    int volume;
    int startpage;
    int endpage;
    Status status;

    public Citation() {
    }

    public Citation(Author[] authors, Status status) {
        this(authors, null, -1, null, -1, -1, -1, status);
        if (status != Status.IN_PREPARATION) {
            throw new CitationException("Only citations in preparation may not contain titles or journals");
        }
    }

    public Citation(Author[] authors, String title, String journal,
                   Status status) {
        this(authors, title, -1, journal, -1, -1, -1, status);
        if (status == Status.PUBLISHED) {
            throw new CitationException("Published citations must have years, volumes and pages");
        }
    }

    public Citation(Author[] authors, String title, int year, String journal, int volume, int startpage, int endpage,
                   Status status) {
        this.authors = authors;
        this.title = title;
        this.year = year;
        this.journal = journal;
        this.volume = volume;
        this.startpage = startpage;
        this.endpage = endpage;
        this.status = status;
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
            builder.append(volume);
            builder.append(", ");
            builder.append(startpage);
            if (endpage > 0) builder.append("-").append(endpage);
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
        builder.append(" <b>").append(volume).append("</b>:");
        builder.append(startpage);
        if (endpage > 0) builder.append("-").append(endpage);
        builder.append("</html>");

        return builder.toString();
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

    class CitationException extends RuntimeException {

        CitationException(String message) {
            super(message);
        }
    }
}
