/*
 * Created on 8 janv. 2007
 */

/** Matrix substitution used for alignment : 
 * transition/transversion matrix for Nucleic Acid ; 
 * identity matrix or Blosum50 matrix for Amino-Acids. 
 * To know if a sequence is Proteic or Nucleic see in LifeSequence. 
 * @see LifeSequence#getSeqType
 * @author Gabriel Chandesris
 */
public class MatrixSub {
    /** The matrix table use in alignment algorithm's*/
    private int[][] matrix;
    /** Determine if the matrix table is "complete" or not (symetric) 
     * @see MatrixSub#aablosum50
     * @see MatrixSub#getValue */
    private boolean complete;
    /** Length of matrix / table ; (ataille = btaille) as it is same number of lines and columns */
    private int ataille,btaille;
    
    /** Default constructor (matrix table is identity matrix for Nucleic Acids) */
    public MatrixSub() {
        this.matrix = this.anidentity();
    }
    /** Indicative Constructor  :
     * <ul><li>0 for identity matrix nucleic acid ; </li>
     * <li>1 for transition-transversion matrix nucleic acid ; </li>
     * <li>2 for blast matrix nucleic acid ; </li>
     * <li>5 for identity matrix for amino-acids ; </li>
     * <li>6 for blosum50 matrix for amino-acids ; </li>
     * <li>7 for pam250 matrix for amino-acids ; </li>
     * <li>...</li></ul>
     * @param a determine the matrix to use */
    public MatrixSub(int a) {
        switch (a) {
        case(0) : this.matrix = this.anidentity();break;
        case(1) : this.matrix = this.antranstrans();break;
        case(2) : this.matrix = this.anblast();break;
        case(5) : this.matrix = this.aaidentity();break;
        case(6) : this.matrix = this.aablosum50();break;
        case(7) : this.matrix = this.aapam250();break;
        }
    }
    
    /** Gives the identity matrix table for Nucleic Acids */
    private int[][] anidentity () {
        int anidentitymatrix[][] = {
            		/* a			t		c		g */
        /* a */ {   1		,	0	,	0	,	0 },
        /* t */ {   0		,	1	,	0	,	0 },
        /* c */ {   0		,	0	,	1	,	0 },
        /* g */ {   0		,	0	,	0	,	1 } };
        complete = true;ataille = btaille = 3;
        return anidentitymatrix;
    }
    /** Gives the transition-transversion matrix table for Nucleic Acid */
    private int[][] antranstrans () {
        int antranstransmatrix[][] = {
        		/* a			t		c		g */
    /* a */ {   3	,	0	,	1	,	0 },
    /* t */ {   0	,	3	,	0	,	1 },
    /* c */ {   1	,	0	,	3	,	0 },
    /* g */ {   0	,	1	,	0	,	3 } };
    complete = true;ataille = btaille = 3;
    return antranstransmatrix;
    }   
    /** Gives the BLAST matrix table (use by blast if in program */
    private int[][] anblast () {
        int anblastmatrix[][] = {
        		/* a			t		c		g */
    /* a */ {   1	,	-3	,	-3	,	-3 },
    /* t */ {   -3	,	1	,	-3	,	-3 },
    /* c */ {   -3	,	-3	,	1	,	-3 },
    /* g */ {   -3	,	-3	,	-3	,	1  } };
    complete = true;ataille = btaille = 3;
    return anblastmatrix;
    }    
    /** Gives the identity matrix table for Amino Acids */
    private int[][] aaidentity() {
	    int aaidentitymatrix[][] = {
		        /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */ 
		/* A */ {  1, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* R */ {  0, 1,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* N */ {  0, 0,	1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* D */ {  0, 0,	0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* C */ {  0, 0,	0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* Q */ {  0, 0,	0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* E */ {  0, 0,	0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* G */ {  0, 0,	0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* H */ {  0, 0,	0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* I */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* L */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* K */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
		/* M */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
		/* F */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
		/* P */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
		/* S */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
		/* T */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
		/* W */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
		/* Y */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 },
		/* V */ {  0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 } 
		        /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */
	};
	    complete = true;ataille = btaille = 19;
	    return aaidentitymatrix;
    }
    /**  Gives the Blosum50 matrix table for Amino Acids ; 
     * ataille and btaille are useful here, but put into identity matrix
     * to make sense and model of methods in this class */
    private int[][] aablosum50() {
	    int aablosum50matrix[][] = {
	        /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */ 
	/* A */ {  5                                                          },
	/* R */ { -2, 7                                                       },
	/* N */ { -1,-1, 7                                                    },
	/* D */ { -2,-2, 2, 8                                                 },
	/* C */ { -1,-4,-2,-4,13                                              },
	/* Q */ { -1, 1, 0, 0,-3, 7                                           },
	/* E */ { -1, 0, 0, 2,-3, 2, 6                                        },
	/* G */ {  0,-3, 0,-1,-3,-2,-3, 8                                     },
	/* H */ { -2, 0, 1,-1,-3, 1, 0,-2,10                                  },
	/* I */ { -1,-4,-3,-4,-2,-3,-4,-4,-4, 5                               },
	/* L */ { -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5                            },
	/* K */ { -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6                         },
	/* M */ { -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7                      },
	/* F */ { -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8                   },
	/* P */ { -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10                },
	/* S */ {  1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5             },
	/* T */ {  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5          },
	/* W */ { -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15       },
	/* Y */ { -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8    },
	/* V */ {  0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5 } 
	        /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */
	    };
	    complete = false;ataille = btaille = 19;
	    return aablosum50matrix;
    }
    /**  Gives the Pam250 matrix table for Amino Acids ; 
     * as for aablosum250 method, attributes are useful here 
     * pam250 is useless because blosum50 exist ! 
     * @deprecated Matrix table pam250 is useless because blosum50 exist ! 
     * @see MatrixSub#aablosum50 */
    private int[][] aapam250() {
	    int aapam250matrix[][] = {
		        /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */ 
		/* A */ {  2                                                          },
		/* R */ { -2, 6                                                       },
		/* N */ {  0, 0, 2                                                    },
		/* D */ {  0,-1, 2, 4                                                 },
		/* C */ { -2,-4,-4,-5, 4                                              },
		/* Q */ {  0, 1, 1, 2,-5, 4                                           },
		/* E */ {  0,-1, 1, 3,-5, 2, 4                                        },
		/* G */ {  1,-3, 0, 1,-3,-1, 0, 5                                     },
		/* H */ { -1, 2, 2, 1,-3, 3, 1,-2, 6                                  },
		/* I */ { -1,-2,-2,-2,-2,-2,-2,-3,-2, 5                               },
		/* L */ { -2,-3,-3,-4,-6,-2,-3,-4,-2, 2, 6                            },
		/* K */ { -1, 3, 1, 0,-5, 1, 0,-2, 0,-2,-3, 5                         },
		/* M */ { -1, 0,-2,-3,-5,-1,-2,-3,-2, 2, 4, 0, 6                      },
		/* F */ { -4,-4,-4,-6,-4,-5,-5,-5,-2, 1, 2,-5, 0, 9                   },
		/* P */ {  1, 0,-1,-1,-3, 0,-1,-1, 0,-2,-3,-1,-2,-5, 6                },
		/* S */ {  1, 0, 1, 0, 0,-1, 0, 1,-1,-1,-3, 0,-2,-3, 1, 3             },
		/* T */ {  1,-1, 0, 0,-2,-1, 0, 0,-1, 0,-2, 0,-1,-2, 0, 1, 3          },
		/* W */ { -6, 2,-4,-7,-8,-5,-7,-7,-3,-5,-2,-3,-4, 0,-6,-2,-5,17       },
		/* Y */ { -3,-4,-2,-4, 0,-4,-4,-5, 0,-1,-1,-4,-2, 7,-5,-3,-3, 0,10    },
		/* V */ {  0,-2,-2,-2,-2,-2,-2,-1,-2, 4, 2,-2, 2,-1,-1,-1, 0,-6, 2, 4 } 
		        /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */
		    };
		    complete = false;ataille = btaille = 19;
		    return aapam250matrix;
    }
    
    /** This class is useful for calling matrix / table in alignment algorithm
     * @param a x-coordinate in table
     * @param b y-coordinate in table
     * @return Value for those (x,y)-coordinate's*/
    public int getValue (int a,int b) {
        int c = 0;int x,y;
        if ( (complete) || (a == b) ) { c = matrix[a][b]; }
        else { 
            if (b >= matrix[a].length) { x = ataille - a;y = btaille - b; }
            else { x = a;y = b; }
            c = matrix[x][y]; 
        }
        return c;
    }
}
