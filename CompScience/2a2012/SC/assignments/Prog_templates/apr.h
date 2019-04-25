/* Structure parm describes the initial string configuration and the
 * simulation/computation parameters.  This structure is filled by
 * the libapr function parseargs or getparms and should not by changed
 * in any way.
 */
extern struct parm {
    enum { NONE=0, SINE=1, PLUCKED=2 } method; /* NONE is an illegal value  */
    int periods;     /* number of sine waves on string in case of SINE      */
    double location; /* global location (0.0 .. 1.0) of plucked position    */
    double stime;    /* time in seconds to simulate                         */
    double delta;    /* time to advance per string computation/iteration    */
    int ntotal;      /* number of points on string in total (global string) */
    int freq;        /* visualize every <freq> iterations, <=0 if none      */
} parms;


extern FILE *visualize_init(int rank);
extern void  visualize (FILE *fp, int  nvalues, double *values, float time);

extern int parseargs(int argc, char *argv[]); /* should be used by INF/AI */
extern int getparms(int argc, char *argv[]);  /* only for BIS students    */

extern int imbalance(MPI_Comm comm, int *npoints, int *nvalues, int *offset);

extern double fit_line(int ndata, double *x, double *y, double *a, double *b);


#ifndef EXIT_FAILURE
# define EXIT_FAILURE 1
# define EXIT_SUCCESS 0
#endif
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#define	TRUE  1
#define FALSE 0
