package dr.app.beauti.util;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.CharSetAlignment;

import java.util.ArrayList;
import java.util.List;

public class CharSet {

    String name;
    List<CharSetBlock> blocks;

    public CharSet(String name) {
        this.name = name;
        blocks = new ArrayList<CharSetBlock>();
    }

    public List<CharSetBlock> getBlocks() {
        return blocks;
    }

    public String getName() {
        return name;
    }

    public void addCharSetBlock(CharSetBlock b) {
        blocks.add(b);
    }

    public Alignment constructCharSetAlignment(Alignment alignment) {

        return new CharSetAlignment(this, alignment);
    }
}
