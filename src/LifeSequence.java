/*
 * Created on 9 janv. 2007
 */
/** This class is used to have an instance for each sequence to put
 * into the phylogenetic tree. 
 * @author Gabriel Chandesris
 */
public class LifeSequence {
    /** Name of the sequence */
    private String name;
    /** The sequence itself */
    private String sequence;
    /** List of abbreviations of Amino-Acids */
    final static public String aminoacids = "ARNDCQEGHILKMFPSTWYV"; 
    /** List of abbreviations of Nitrogen Bases (ADN) */
    final static public String nitrogenbases = "atcg";
    
    /** This constructor make the easiest work 
     * @param nam For the name of the sequence
     * @param seq For the sequence*/
    public LifeSequence (String nam,String seq) {
        this.name = nam;this.sequence = seq;
    }
    
    /** Return the name of the sequence */
    public String getName() { return name; }
    /** return the sequence */
    public String getSequence() { return sequence; }
    
    /**
     * Determine if sequence is Nucleic Acids sequence (i.e. ADN, 
     * but not ARN because of unpresent uracil)
     *  or a Amino-Acids sequence (i.e. a protein). 
     * @param seqtest
     * @return false if Nucleic ; true if Proteic
     */
    public static boolean getSeqType (String seqtest) {
        boolean isproteic = false;

        for (int i = 0 ; i < 20 ; i++) {
            if (seqtest.charAt(0) == aminoacids.charAt(i)) {
                isproteic = true; } }
        // optionnal part is following
        // because if it is not proteic (true), it is nucleic (false)
        for (int i = 0 ; i < 4 ; i++) {
            if (seqtest.charAt(0) == nitrogenbases.charAt(i)) {
                isproteic = false; } }
        return isproteic;
    }

    /** Return length of the sequence itself */
    public int getLength() { return this.sequence.length(); }
    
    /** To get the String of name and sequence. */
    public String getAff() {
        return this.name+"\n"+this.sequence+"\n";
    }
    
    /** To get only the name on String.  */
    public String toString() {
        String txt_aff = "";
        if (this.name.length() > 25) { txt_aff = this.name.substring(0,25); }
        else { txt_aff = this.name; }
        return txt_aff+"\n"; 
    }
}
