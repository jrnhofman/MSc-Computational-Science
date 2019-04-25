#ifndef MYAPR_H
#define MYAPR_H 

/* Structure parm describes the initial string configuration and the
 * simulation/computation parameters.  This structure is filled by
 * the libapr function parseargs or getparms and should not by changed
 * in any way.
 */
extern struct parm {
    double threshold;    /* time in seconds to simulate                         */
    double delta;    /* time to advance per string computation/iteration    */
    int ntotal;      /* number of points on string in total (global string) */
    int freq;        /* visualize every <freq> iterations, <=0 if none      */
} parms;

//Parsing functions
extern int parseargs(int argc, char *argv[]); /* should be used by INF/AI */

//Get command line arguments; rank 0 reads them from stdin, and then broadcasts
// the values to all other processes. This function makes use of parseargs.
extern int getparms(int argc, char *argv[]);  /* only for BIS students    */

// NOTE: below is unused I think
#ifndef EXIT_FAILURE
# define EXIT_FAILURE 1
# define EXIT_SUCCESS 0
#endif

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#define	TRUE  1
#define FALSE 0

#endif
