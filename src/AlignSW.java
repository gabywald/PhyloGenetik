/*
 * Created on 8 janv. 2007
 */
/** In this class, application of the Smith and Waterman algorithm
 * Optimization of Needelman and Wunsch Algorithm
 * @author Gabriel Chandesris
 */
public class AlignSW extends Alignment {
    /** Special attribute of this class for maximum local value (3 attributes) */
    private int maxvalue = 0,imaxvalue = 0,jmaxvalue = 0;
    
    /** Constructor of object (no default constructor) is part 1 of the SW algorithm 
     * in the objective to get the score. Part 2 is on method alignment. Default gap is 6.
     * @param squn first sequence to align
     * @param sqde second sequence to align
     * @see AlignSW#alignment
     * @see AlignNW
     */
    public AlignSW (LifeSequence squn,LifeSequence sqde)
		throws DiffTypSeqExc {
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
    /** Constructor of object (no default constructor) is part 1 of the SW algorithm 
     * in the objective to get the score. Part 2 is on method alignment. 
     * Similar to constructor of NW constructor
     * @param squn first sequence to align
     * @param sqde second sequence to align
     * @param gap cost of gap (-)
     * @see AlignSW#alignment
     * @see AlignNW
     */
    public AlignSW (LifeSequence squn,LifeSequence sqde,int gap)
		throws DiffTypSeqExc {
        this.setSequn(squn);this.setSeqde(sqde);this.setGap(gap);
        this.matrixscore = new int[squn.getLength()+1][sqde.getLength()+1];
        this.alignment = new boolean[squn.getLength()+1][sqde.getLength()+1];
        
        try { this.buildclass(squn,sqde,gap); }
        catch (DiffTypSeqExc e) { this.buildexcep(-1); }
        

    }
    /** This constructor is calling when error appear. 
     * Initialize the matrix of alignment to a default value (-1)
     * @param a The value to give ; -1 is a good solution ??
     */
    public AlignSW (int a) { this.buildexcep(a); }
    
    /** Smith and Waterman alignment Part 2 ;
     * TODO 
     * @deprecated This method has not be done !
     * @see Alignment#algorithm
     */
    public void algorithm() {
        int i = this.getSequn().length();
        int j = this.getSeqde().length();
        
        this.alignment[this.getSeqde().length()]
                       [this.getSequn().length()]                          
                        = false;
        // false until it's not done !
    }
    
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
        for (int i = 0 ; i < squn.getLength() ; i++) {
            this.matrixscore[0][i] = 0;
            this.alignment[0][i] = false;
            for (int j = 0 ; j < squn.getLength() ; j++) {
                this.matrixscore[j][0] = 0;
                this.alignment[j][0] = false;this.alignment[j][i] = false;
                this.matrixscore[j][i] = max (0,
                        this.matrixscore[j-1][i-1]-this.calculscore(j,i),
                        this.matrixscore[j-1][i]-this.getGap(),
                        this.matrixscore[j][i-1]-this.getGap());
                if (this.matrixscore[j][i] > maxvalue) {
                    this.maxvalue = this.matrixscore[j][i];
                    this.imaxvalue = i;this.jmaxvalue = j;
                }
            }
        }
        this.setScore(this.maxvalue);
    }
}
