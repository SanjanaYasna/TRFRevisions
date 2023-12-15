/*Minimalistic, no reporting results, at alll!*/
#ifndef TRFRUN_NEW
#define TRFRUN_NEW
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>
#include "tr30dat.h"
#include "tr30dat.c"
#include "trfclean.h"


/* These declarations moved by Yevgeniy Gelfand on Jan 27, 2010  */
/* To have smaller sequences not send results */
/* to disc to improve performance             */

/* define Index List structure */

struct index_list
{
	int     count;      /* indicates order in original file */
	char    ref[45];    /* records label for linking */
	int     first;      /* first index */
	int     last;       /* last index */
	int     period;     /* period size */
	float   copies;     /* number of copies */
	int     size;  /* consensus size */
	int     matches;
	int     indels;
	int     score;
	int     acount;
	int     ccount;
	int     gcount;
	int     tcount;
	float   entropy;
	char*   pattern;

	struct index_list* next;
};
typedef struct index_list IL;

IL *GlobalIndexList = NULL;
IL *GlobalIndexListTail = NULL;

void    FreeList(IL * headptr);

#ifndef _MAX_PATH
#define _MAX_PATH 260
#endif

// how much memory to allocate at first when loading a sequence
// (added by Eugene Scherba, 2010-02-16)
#define MEMBLOCK (10*1024*1024)


//function definitions
int LoadSequenceFromFileEugene(FASTASEQUENCE *pseq,FILE* fp); /* may use stdin and better file reading for handling large files over 2GB */
void TRFControlRoutine(void);
void TRF(FASTASEQUENCE* pseq);
void PrintError(char* errortext);
void PrintProgress(char* progresstext);
void SetProgressBar(void);


#ifdef UNIXGUI
extern void set_progress_bar( double fraction );
#endif

/*******************************************************
 *   loads a sequence from an open file. If the routine
 *   comes to another sequence header while reading a
 *   sequence the sequence returns 1 to indicate that
 *   more sequences are present. Otherwise the routine
 *   returns 0 to indicate EOF or -1 to indicate error.
 *   The member sequence must be NULL before the routine
 *   iscalled. The calling function must free the allo-
 *   cated memory after use.
 ********************************************************/

int LoadSequenceFromFileEugene(FASTASEQUENCE *pseq, FILE* fp)
{
	int i, j, c;
	int next = -1;	// whether a next sequence was encountered
	char *ptemp;
	char to_upper;

	// read the FASTA '>' symbol
	c = getc(fp);
	if ((char)c != '>' || c == EOF) return -1; // invalid format

	// read name and description text
	for(i = 0; i < MAXSEQNAMELEN - 1; i++) {
		c = getc(fp);
		if (c == 10 || c == 13) {
			break;
		} else if (c == EOF) {
			PrintError("FASTA input terminated too early");
			return -1;
		} else {
			pseq->name[i] = (char)c;
		}
	}
	pseq->name[i] = '\0';

	// if line was not read completely flush the rest
	if (i == MAXSEQNAMELEN - 1) {
		c = 0;
		while (c != 13 && c != 10 && c != EOF) c = getc(fp);
	}

	// read sequence into buffer
	pseq->sequence = NULL;
	for (i = 0; i < 26; i++) pseq->composition[i] = 0;
	to_upper = 'A' - 'a';

	i = 0;
	for (j = 0; c != EOF && (char)c != '>'; j += MEMBLOCK) {
		if ((ptemp = realloc(pseq->sequence, sizeof(char)*(i + MEMBLOCK + 1))) == NULL) {
			PrintError("Insufficient memory");
			return -1;
		} 
		pseq->sequence = ptemp;
		for(i = j; i < j + MEMBLOCK;) {
			c = getc(fp);			// get a character from file
			if ((char)c >= 'A' && (char)c <= 'Z') { // in upper-case range of alpha characters
				pseq->sequence[i] = (char)c;
				pseq->composition[c - (int)'A']++;
				i++;
			} else if ((char)c >= 'a' && (char)c <= 'z') {	// in lower-case range of alpha characters
				pseq->sequence[i] = (char)c + to_upper;	// make upper-case
				pseq->composition[c - (int)'a']++;
				i++;
			} else if ((char)c == '>') {	// break if another sequence found
				next = 1;
				ungetc('>', fp);
				break;
			} else if (c == EOF) {		// break if end of file
				next = 0;
				break;
			}
		}
	}
	pseq->length = i;		// set sequence length
	if (i > 0) {
		pseq->sequence[i] = '\0';	// terminate sequence text as a string
    }
	// compute member
	pseq->nucleotides =
		pseq->composition['A'-'A'] + pseq->composition['C'-'A'] +
		pseq->composition['G'-'A'] + pseq->composition['T'-'A'];
	return next;
}



#endif 
