/*
 * Created on 4 janv. 2007
 */

/** 
 * This abstract class is use to define correctly attributes and methods
 * of alignments like NW and SW (distinct classes which extends this one). 
 * See the Direct subclasses for more details. 
 * @author Gabriel Chandesris
 */
public abstract class Alignment {
    /** Table made to build score. */
    protected int[][] matrixscore;
    /** Table of alignment. */
    protected boolean[][] alignment;
    /** Final score of alignment. */
    protected int score;
    /** Value of gap. */
    protected int gap;
    /** Value of deafult gap (simple constructor) */
    protected static int gap_default = 6;
    /** Sequence. */
    protected LifeSequence seq1, seq2;
    /* DONE substituion matrix
     * tansition/transversion for nucleic acid
     * blosum50 for amino-acids / protein*/
    /** Substitution matrix table */
    protected MatrixSub substitutiontable;
    /** Protein alignment or Nucleic Alignment. */
    protected boolean isprotein;
    /**  */
    
    /** Default constructor to set the proteic ans substitution table.  */
    public Alignment () { ; }

    
    /** Algorithm of alignment must be defined in inheritant classes
     * with action on coutalignment[][] with attributes 
     * score, gap, seq1 and seq2 */
    public abstract void algorithm ();
    
    /** This methode is used by the constructors of this class */
    protected abstract void buildclass 
    		(LifeSequence squn,LifeSequence sqde,int gap) 
    		throws DiffTypSeqExc;
    
    /** If an exception appear in execution of the constructor of this class. */
    protected void buildexcep (int a) {
        this.setSequn(null);this.setSeqde(null);this.setGap(0);
        // create the instance of matrix score
        this.matrixscore = new int[this.getSequn().length()][];
        for (int i = 0 ; i < this.matrixscore.length ; i++) {
            this.matrixscore[i] =  new int[this.getSeqde().length()]; }
        // create the instance of alignment
        this.alignment = new boolean[this.getSequn().length()][];
        for (int i = 0 ; i < this.alignment.length ; i++) {
            this.alignment[i] = new boolean[this.getSeqde().length()]; }
        for (int i = 0 ; i < this.getSequn().length() ; i++) {
            for (int j = 0 ; j < this.seq2.getLength() ; j++) {
                this.matrixscore[i][j] = a;
                this.alignment[i][j] = false; } }
    }
    
    /** This method build the alignment and return it
     * @deprecated Method not done : don't need it here for now
     * @return A three-lines element to be show */
    public String getAlignment() {
        String line1="",line2="",line3="";
        if (alignment[seq1.getSequence().length()][seq2.getSequence().length()]) {
            //TODO catch true values
        }
        else {
            //TODO exception because alignment not done
        }
        return line1+"\n"+line2+"\n"+line3+"\n";
    }
    
    /** This method is here only to get the value in object MatrixSub
     * @see MatrixSub#getValue  */
    public int getSubValue(int a,int b) { 
        return this.substitutiontable.getValue(a,b);
    }
    
    /**This method is used for comparison between two char in the two sequences
     * @param j is reference to a char of the first sequence
     * @param i is reference to a char of the second sequence
     * @return 1 for match and 0 for mismatch
     * @deprecated Method calculscore is a better way to do and give a result between two chars...
     * @see Alignment#calculscore
     */
    public int calculmatchmismatch(int j,int i) {
        int match = 0;
        if (seq1.getSequence().charAt(j) == seq2.getSequence().charAt(i)) { match = 1; }
        return match;
    }
    
    /** This method is a better version of calcumatchmismatch, and refer 
     * to some elements of LifeSequence to compute the score. 
     * @see LifeSequence#aminoacids
     * @see LifeSequence#nitrogenbases
     * @param i is reference to a char of the first sequence
     * @param j is reference to a char of the second sequence
     * @return Score in the substitution matrix table for those two char's
     */
    public int calculscore(int i,int j) {
        int a = 0;int b = 0;int k = 0;String charSequences;
        if (isprotein) { charSequences = LifeSequence.aminoacids;k = 19; } 
        else { charSequences = LifeSequence.nitrogenbases;k = 3; }
        // k = charSequences.length(); // 20 or 4
        while ( (a<k) && (seq1.getSequence().charAt(i) != charSequences.charAt(a))) { a++; }
        while ( (b<k) && (seq2.getSequence().charAt(j) != charSequences.charAt(b))) { b++; }
        return this.getSubValue(a,b);
    }
    
    /** Return score */
    public int getScore() { return this.score; }
    /** Return gap */
    public int getGap() { return this.gap; }
    /** Return first sequence (string) */
    public String getSequn(){ return this.seq1.getSequence(); }
    /** Return second sequence (string) */
    public String getSeqde(){ return this.seq2.getSequence(); }
    /** Return length of first sequence. */
    public int getSequnLength() { 
        if (this.seq1 == null) { return 0; }
        else { return this.getSequn().length(); }
    }
    /** Return length of second sequence. */
    public int getSeqdeLength() { 
        if (this.seq2 == null) { return 0; }
        else { return this.getSeqde().length(); }
    }
    
    /** Return String to input into the textarea of interface. */
    public String getAff () {
        return this.seq1.getAff()+this.seq2.getAff()+"Gap : "+gap+"\nScore : "+score+"\nScore Matrix : "+this.getAffMatrixScores();
    }
    
    /** Return the matrix of score to be put into the textarea of interface.  */
    public String getAffMatrixScores() { 
        String aff_matrixscore = "";
        String linebreak = "";
        int k = 0;
        while (k < this.getSequnLength()) { 
            linebreak = linebreak+"----";
            k++; 
        }
        aff_matrixscore = aff_matrixscore+"\n"+linebreak+"\n |";
        for (int i = 0 ; i < this.getSequnLength() ; i++) {
            for (int j = 0 ; j < this.getSeqdeLength() ; j++) {
                String temp = matrixscore[i][j]+"";
                String insert = "|";
                if (temp.length() == 3) { insert = temp+insert; }
                if (temp.length() == 2) { insert = " "+temp+insert; }
                if (temp.length() == 1) { insert = " "+temp+" "+insert; }
                aff_matrixscore = aff_matrixscore+insert;
            }
            aff_matrixscore = aff_matrixscore+"\n"+linebreak+"\n |";
        }
        return aff_matrixscore; }
    
    /** Return a String value for this object. */
    public String toString() {
        String aff_txt="";
        aff_txt = aff_txt+"Gap : "+this.gap+"\n";
        aff_txt = aff_txt+"Type : "+((this.isprotein)?"proteic":"nucleic")+"\n";
        aff_txt = aff_txt+"First sequence : \n"+this.seq1.toString()+"\n";
        aff_txt = aff_txt+"Second sequence : \n"+this.seq1.toString()+"\n";
        aff_txt = aff_txt+"Score of Alignment : "+this.score+"\n";
        return aff_txt;
    }

    // /** Return first sequence (LifeSequence) */
    // public LifeSequence getSequn(){ return seq1; }
    // /** Return second sequence (LifeSequence) */
    // public LifeSequence getSeqde(){ return seq2; }
    
    /** Set the first sequence with given argument */
    public void setSequn(LifeSequence squn){ this.seq1 = squn; }
    /** Set the second sequence with given argument */
    public void setSeqde(LifeSequence sqde){ this.seq2 = sqde; }
    /** Set the score  with given argument*/
    public void setScore(int scr) { this.score = scr; }
    /** Set the gap  with given argument*/
    public void setGap (int gp) { this.gap = gp; }
    /** Set the default gap to a new value */
    public static void setGapDefault (int gp) { Alignment.gap_default = gp; }
    /** Minimum of two integers */
	static int min(int a,int b){
	    if (a < b) return a;
	    else return b;
	}
    /** Minimum of three integers */
	static int min(int a, int b, int c) { return min(min(a,b),c); }
	/** Minimum of four integers */
	static int min(int a, int b, int c,int d) { return min(min(a,b),min(c,d)); }
	/** Maximum of two integers */
	static int max(int a,int b){
	    if (a > b) { return a; }
	    else return b;
	}
    /** Maximum of three integers */
	static int max(int a, int b, int c) { return max(max(a,b),c); }
	/** Maximum of four intergers */
	static int max(int a, int b, int c,int d) { return max(max(a,b),max(c,d)); }


	/*
	 * Created on 17 janv. 2007
	 */
	/** This exception appear when the two sequences are Nucleic and Proteic.
	 * Throw in constructor and subclasse's, catch in Subclasses. 
	 * @see Alignment#Alignment
	 * @see AlignNW
	 * @see AlignSW
	 * @author Gabriel Chandesris
	 */
	public class DiffTypSeqExc extends Exception { ; }
}


