/* FASTASPLITN.C
Version:   1.0.4  17-DEC-2007
Author:    David Mathog, Biology Division, Caltech
email:     mathog@caltech.edu
Copyright: 2001 David Mathog and California Institute of Technology (Caltech)

Description:
   Splits a fasta file into a N files, with cycling entries among the N
   output  files. This does a better job of randomizing the split content
   than did FASTASPLIT.C, which doesn't work well when there is a systematic
   variation in the input file.  For instance, in the NCBI "nt" database the
   entries get  longer and longer from the beginning to the end, so that
   regular FASTASPLIT  produced a final file 5 times bigger than the first one
   (all with the same  number of sequences though!)

   If the optional P value is specified, it emits output to
   stdout with contents identical to the Pth fragment of a full
   N file set. 

   If in addition the optional C value is specified the number of sequential
   entries emitted to each stream may be changed from the default of 1.

   Very little error checking.
   
   No input line may exceed 1M  characters.
   
   Compiles cleanly with:
   
   gcc -Wall -ansi -pedantic -o fastasplitn fastasplitn.c

Arguments are:

 1:  infile (- or STDIN mean read from stdin)
 2:  N number of files to produce. 
 3:  P phase. If specified and in range 1-N only a single
       file  is produced as if it was the Pth of the N files. 
       If P=0 then N files are produced.  Anything else is an error.
 4:  C cycle.  If specified must be >=1.  If C=3, N=2 then
      1,2,3,7,8,9, etc.    -> file1
      4,5,6,10,11,12, etc. -> file2
      C=1 is the default.
      
 
License terms:
    You may run this program on any platform. You may
    redistribute the source code of this program subject to
    the condition that you do not first modify it in any way.
    You may  distribute binary versions of this program so long
    as they were compiled from unmodified source code.  There
    is no charge for using this software.  You may not charge
    others for the use of this software.

Bug reports:  Please report bugs via email.

1.0.4  Slightly clarified command line usage, after 
       a suggestion by Bernd Brandt.
1.0.3  Added code to implement C value.  Useful for splitting files
       into sequential sequences if one knows ahead of time how many
       sequences are present in the input file.  Ie, if there are 20
       and N=4 C=5 P=0 then 4 files will be produced with the first
       holding seqs 1->5, the second 6->10, and so forth.
1.0.2  Added code to read symbol "SPLITFRAGTEMPLATE" which replaces
       default value template for output filename creation.  Cleaned
       up printf so that all messages go to stderr.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>  /* for toupper */
#define MYMAXSTRING 100000
#define MAXOUTFILES 1000  /* WAY more than should ever be needed - or possible! */


int  lcl_strcasecmp(char *s1, char *s2){
int c1;
int c2;
  for(; ;s1++,s2++){
    c1=toupper(*s1);
    c2=toupper(*s2);
    if(c1 < c2)return -1;
    if(c1 > c2)return  1;
    if(c1 == 0)return  0;  /*c2 also is 0 in this case */
  }
}

int main(int argc, char *argv[]){
char *infile;
char bigstring[MYMAXSTRING];
char outname[200];
char *thetemplate;
int  n,p,c,count,fragcount,cyccount;
FILE *fin;
FILE *fout[MAXOUTFILES]; 
  
  fragcount=0;
  count=0;
  n=0;
  p=0;
  c=1;
  if( (argc < 3)                                   ||
      (argc > 6)                                   ||
      (argc >=3 && (sscanf(argv[2],"%d",&n)==EOF)) ||
      n < 1                                        ||
      (argc >=4 && 
        (     (sscanf(argv[3],"%d",&p)==EOF)  ||
              p < 0                           ||
              p > n
        )
      )                                            ||
      (argc >=5 && 
        (     (sscanf(argv[4],"%d",&c)==EOF)  ||
              c < 1
        )
      )
    ){
    (void) fprintf(stderr,"Fastasplitn: invalid input parameters: n %d p %d c %d\n\n", n, p, c);
    (void) fprintf(stderr,"Function:  Splits a fasta file into N data streams. \n");
    (void) fprintf(stderr,"           Default is to send each stream to a separate output file.\n");
    (void) fprintf(stderr,"           Optional command line parameters are indicated with [].\n\n");
    (void) fprintf(stderr,"Usage:     fastasplitn infile N [P [C]]\n");
    (void) fprintf(stderr,"  infile = input file.  \"-\" or \"stdin\" = read from stdin\n");
    (void) fprintf(stderr,"  N =      number of data streams to produce, N must be >0\n");
    (void) fprintf(stderr,"  P =      0 or omitted, emit N streams to N files\n");
    (void) fprintf(stderr,"    =      1->N, only emit contents of Pth stream (1-N), send to stdout\n");
    (void) fprintf(stderr,"  C =      emit C sequences in order for each stream, C>=1, 1 is default\n\n");
    (void) fprintf(stderr,"Examples: in all cases input has 20 sequences\n\n");
    (void) fprintf(stderr,"  N=4 P=0 C=5:\n");
    (void) fprintf(stderr,"    1->5 to file1, 6->10 to file2, 11->15 to file3, 16->20 to file4\n");
    (void) fprintf(stderr,"  N=4 P=2 C=3:\n");
    (void) fprintf(stderr,"    4->6,16->17 to stdout\n\n");
    (void) fprintf(stderr,"Default fragment file name template is \"frag%%3.3d\" \n");
    (void) fprintf(stderr,"The string in the symbol SPLITFRAGTEMPLATE overrides this default\n");
    exit(EXIT_FAILURE);
  }
  if((n < 1) || (n > MAXOUTFILES)){
    (void) fprintf(stderr,"The number of files must be between 1 and %d\n",MAXOUTFILES);
    exit(0);
  }

  infile=argv[1];
  
  if(p==0)(void) fprintf(stderr,"Processing %s with n=%d p=%d c=%d\n",infile,n,p,c);
  
  if(lcl_strcasecmp(infile,"-") == 0 || lcl_strcasecmp(infile,"stdin") == 0){
    fin=stdin;
  }
  else {
    fin=fopen(infile,"r");
  }
  if(fin==NULL){
   (void) fprintf(stderr,"Could not open input file %s\n",infile);
    exit(0);
  }

  /* see if template is to use default or from a symbol */
  
  thetemplate=getenv("SPLITFRAGTEMPLATE");
  if(thetemplate==NULL){
     thetemplate=getenv("splitfragtemplate");
     if(thetemplate==NULL){
       thetemplate=malloc(30);
       (void) strcpy(thetemplate,"frag%.3d");
    }
  }
  if(p==0)(void) fprintf(stderr,"the template is %s\n",thetemplate);


  /* open all the output files */

  for(fragcount=0; fragcount<n ;fragcount++){
     (void) sprintf(outname,thetemplate,fragcount+1);
     if(p==0){
       (void) fprintf(stderr,"Opening output file %s\n",outname);
       fout[fragcount] = fopen(outname,"w");
       if(fout[fragcount]==NULL){
          (void) fprintf(stderr,"Could not open output file %s\n",outname);
          exit(EXIT_FAILURE);
       }
     }
     else if(p==fragcount+1){
       fout[fragcount] = stdout;
     }
  }

  /* now send fasta entries to each one */

  fragcount=0;
  cyccount=0;
  count = 0;
  while( fgets(bigstring,MYMAXSTRING,fin) != NULL){
    if(bigstring[0] == '>'){
      count++;
      cyccount++;
      if(cyccount > c){
        cyccount=1;
        fragcount++;
        if(fragcount >= n)fragcount=0;
      }
    }
    if(p==0 || p==fragcount+1)(void) fprintf(fout[fragcount],"%s",bigstring);
  }
  if(p==0)(void) fprintf(stderr,"All done, entries processed: %d\n",count);
  exit(EXIT_SUCCESS);
}
