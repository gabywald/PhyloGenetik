/*
 * Created on 4 janv. 2007
 */

/**
 * In this class, application of the Needelman and Wunsch algorithm
 * @author Gabriel Chandesris
 */
public class AlignNW extends Alignment {

    /**
     * Constructor of object (no default constructor) is part 1 of the NW algorithm 
     * in the objective to get the score. Part 2 is on method alignment. Default gap is 6.
     * @param squn first sequence to align
     * @param sqde second sequence to align
     * @see AlignNW#alignment
     * @see AlignSW
     */
    public AlignNW (LifeSequence squn,LifeSequence sqde) {
        this.setSequn(squn);this.setSeqde(sqde);this.setGap(gap_default);
        // create the instance of matrix score
        this.matrixscore = new int[this.getSequnLength()][];
        for (int i = 0 ; i < this.matrixscore.length ; i++) {
            this.matrixscore[i] =  new int[this.getSeqdeLength()]; }
        // create the instance of alignment
        this.alignment = new boolean[this.getSequnLength()][];
        for (int i = 0 ; i < this.alignment.length ; i++) {
            this.alignment[i] = new boolean[this.getSeqdeLength()]; }
        
        try { this.buildclass(squn,sqde,gap_default); }
        catch (DiffTypSeqExc e) { this.buildexcep(-1); }
    }
    /**
     * Constructor of object (no default constructor) is part 1 of the NW algorithm 
     * in the objective to get the score. Part 2 is on method alignment.
     * @param squn first sequence to align
     * @param sqde second sequence to align
     * @param gap cost of gap (-)
     * @see AlignNW#alignment
     * @see AlignSW
     */
    public AlignNW (LifeSequence squn,LifeSequence sqde,int gap) {
        this.setSequn(squn);this.setSeqde(sqde);this.setGap(gap);
        // create the instance of matrix score
        this.matrixscore = new int[this.getSequnLength()][];
        for (int i = 0 ; i < this.matrixscore.length ; i++) {
            this.matrixscore[i] =  new int[this.getSeqdeLength()]; }
        // create the instance of alignment
        this.alignment = new boolean[this.getSequnLength()][];
        for (int i = 0 ; i < this.alignment.length ; i++) {
            this.alignment[i] = new boolean[this.getSeqdeLength()]; }
        
        try { this.buildclass(squn,sqde,gap); }
        catch (DiffTypSeqExc e) { this.buildexcep(-1); }
    }
    /** This constructor is calling when error appear. 
     * Initialize the matrix of alignment to a default value (-1)
     * @param a The value to give ; -1 is a good solution ??
     */
    public AlignNW (int a) { this.buildexcep(a); }
    
    /** This methode is used by the constructors of this class */
    protected void buildclass (LifeSequence squn,LifeSequence sqde,int gap) 
    		throws DiffTypSeqExc {
        // initialize this instance of object class NW
        // this.setSequn(squn);this.setSeqde(sqde);this.setGap(gap);
        boolean isprot1 = LifeSequence.getSeqType(this.getSequn());
        boolean isprot2 = LifeSequence.getSeqType(this.getSeqde());
        // DONE exception if the two LifeSequence's are of 
        // both different Nucleic and Proteic types...
        if (isprot1 != isprot2) { throw new DiffTypSeqExc(); }
        else {
            this.isprotein = isprot1;
            if (isprotein) { substitutiontable = new MatrixSub(6); } 
            else { substitutiontable = new MatrixSub(1); }
        }

        // making alignment part 1 : comparison between the two sequences
        // start after first line and column because of initializing values
        for (int i = 1 ; i < this.getSequn().length() ; i++) {
            this.matrixscore[i][0] = i*this.getGap();
            this.alignment[i][0] = false;
            for (int j = 1 ; j < this.seq2.getLength() ; j++) {
                this.matrixscore[0][j] = j*this.getGap();
                this.alignment[0][j] = false;
                this.alignment[i][j] = false;
                this.matrixscore[i][j] = min (
                        this.matrixscore[i-1][j-1]+this.calculscore(i,j),
                        this.matrixscore[i-1][j]+this.getGap(),
                        this.matrixscore[i][j-1]+this.getGap());
            }
        }
        this.setScore(this.matrixscore[this.getSequn().length()-1][this.getSeqde().length()-1]);
    }
    
    /**Needelman and Wunsch alignment Part 2
     * @see Alignment#algorithm */
    public void algorithm() {
        int i = this.seq1.getLength();
        int j = this.seq2.getLength();
        int a,b,c;
        while ((i!=0)&&(j!=0)) {
            if (this.getSequn().charAt(i) == this.getSeqde().charAt(j)) {
                i--;j--;this.alignment[j][i] = true;
            }
            else {
                a = this.matrixscore[j-1][i-1];
                b = this.matrixscore[j][i-1];
                c = this.matrixscore[j-1][i];
                if (max(a,b,c) == a) { j--;i--; }
                if (max(a,b,c) == b) { i--; }
                if (max(a,b,c) == c) { j--; }
                this.alignment[j][i] = true;
            }
        }
        this.alignment[this.getSeqde().length()][this.getSequn().length()] = true;
    }

}
