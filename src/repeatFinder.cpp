#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "repeatFinder.h"

using namespace std;

#define REPEAT_LEN                    11
#define HASH_LEN     (1<<(REPEAT_LEN*2))

unsigned INIT_DNA_LEN = 120000000; 

char inputfile[256] = "test.fasta";
char *dna;
int converter[128], complement[128];
int dna_len;
int output_rep_len = REPEAT_LEN;
int debug = 0;

vector<int> allrepeats[HASH_LEN];
struct repeat{
	int fst;
	int sec;
	int len;
	int seclen;
	int visited;
	int exact;
}R;

vector<repeat> rep;
int gap_len = 0; //allow gap < gap_len
int ppno = 0; // which prophage are we working on

void input()
{
	FILE *f;
	dna = new char[INIT_DNA_LEN+1];
	if(!dna){
		printf("Can not allocate memory\n");
		exit(0);
	}

	//input dna
	f = fopen(inputfile, "r");
	if(f == NULL){
		printf("Can not find input file %s\n",inputfile);
		exit(0);
	}
	char *my_new_check = fgets(dna,INIT_DNA_LEN,f);
	(void) my_new_check; // remove an unused parameter warning
	
	strcpy(dna,"");
	while(fscanf(f, "%s", dna+dna_len)==1) {
		int t = strlen(dna+dna_len);
		dna_len += t;
	}
	fclose(f);

	//printf("dna len = %d\n",dna_len);
}

// find all the 11 length substring and store their start position
void find_repeats()
{
	int key,start,keylen=2*(REPEAT_LEN-1);

	key = 0;
	for(start = 0;start < REPEAT_LEN; start++)
		key = (key<<2) + converter[(int) dna[start]];
	allrepeats[key].push_back(0);

	for(start = 1;start < dna_len-REPEAT_LEN+1; start++)
	{
		key = ((key&((1<<keylen)-1))<< 2) + converter[(int) dna[start+REPEAT_LEN-1]];
		allrepeats[key].push_back(start);
	}

	//find reverse repeat
	key = 0;
	for(start = dna_len-1;start >dna_len-1-REPEAT_LEN; start--)
		key = (key<<2) + converter[complement[(int) dna[start]]];
	allrepeats[key].push_back((dna_len-1)*(-1));
	
	for(start = dna_len-2;start >REPEAT_LEN-2; start--)
	{
		key= ((key&((1<<keylen)-1))<<2) + converter[complement[(int) dna[start-REPEAT_LEN+1]]];
		allrepeats[key].push_back(start*(-1));
	}
}

void find_maxlen(int fst, int sec)
{
	int i,j,k;
	//  check whther to compute //
	
	if(sec-fst < REPEAT_LEN)
		return;
	
	if (fst >0)
		if(dna[fst-1] == dna[sec-1])
			return;

	k=0;
	for(i = REPEAT_LEN+fst, j = REPEAT_LEN+sec ; i<sec && j<dna_len; i++, j++)
		if (dna[i] == dna[j])
			k++;
		else 
			break;

    if (debug && sec > dna_len) {
        fprintf(stderr, "In find_maxlen we added a repeat with sec %d that is longer than dna_len (%d)\n", sec, dna_len);
    }

	R.fst=fst+1;
	R.sec=sec+1;
	R.len=REPEAT_LEN + k;
	R.seclen = R.len;
	R.visited = 0;
	R.exact = 0;
	rep.push_back(R);
}

void find_maxlen_rev(int fst, int sec)
{
	int i,j,k;
	
	sec = sec *(-1);
	//  check whther to compute //
	if(sec-fst +1< 2 * REPEAT_LEN)
		return;
	
	if (fst >0 && sec < dna_len-1 )
		if(dna[fst-1] == complement[(int) dna[sec+1]])
			return;
///////////////

	k=0;
	for(i = REPEAT_LEN+fst, j = sec - REPEAT_LEN ; i<dna_len && j>-1 && i<j; i++, j--)
		if (dna[i] == complement[(int) dna[j]])
			k++;
		else 
			break;
	
	R.fst=fst+1;
	R.sec=sec*(-1)-1;
	R.len=REPEAT_LEN + k;
	R.seclen = R.len;
	R.visited = 0;
	R.exact = 0;
	rep.push_back(R);
}

void extend_repeats()
{
	int i,j,key,keylen=2*(REPEAT_LEN-1);

	key = 0;
	for(i = 0;i < REPEAT_LEN; i++)
		key = (key<<2) + converter[(int) dna[i]];
	for(j=0; j < (int) allrepeats[key].size(); j++)
		if(allrepeats[key][j]<0)
			find_maxlen_rev(0,allrepeats[key][j]);
		else
			find_maxlen(0,allrepeats[key][j]);

	for(i =1;i<dna_len-REPEAT_LEN+1;i++){
		key = ((key&((1<<keylen)-1))<< 2) + converter[(int) dna[i+REPEAT_LEN-1]];
		for(j=0; j < (int) allrepeats[key].size(); j++)
			if(allrepeats[key][j]<0)
				find_maxlen_rev(i,allrepeats[key][j]);
			else
				find_maxlen(i,allrepeats[key][j]);
	}
}

int check_extend(int fst,int n)
{
	int i,j = rep[fst].fst + rep[fst].len -1 +n, len = rep.size(),k, head , tail, mid;


	// binery search
	head = fst+1;
	tail = len-1;
	mid = (head +tail)/2;
	while(rep[mid].fst!= j && head<=tail){
		mid = (head +tail)/2;
		if(rep[mid].fst<j)
			head = mid+1;
		else
			tail = mid-1;
	}

	if (head <= tail)
		return -1;

	i = mid-1;
	// while (rep[mid].fst == rep[i].fst && i >-1)
	while(rep[mid].fst == rep[i].fst && i > tail-1)
		i--;
	i++;
////
	k = rep[fst].sec+rep[fst].seclen -1;
	
	for(; i < (int) rep.size() && rep[i].fst == j;i++)
		if(rep[i].visited == 0){
			if ( (rep[i].fst + rep[i].len -1) <rep[fst].sec || (rep[i].fst + rep[i].len -1)<(rep[i].sec)*(-1)-rep[i].seclen+1)
				//check for 2nd copy
				if(rep[i].sec-k <= gap_len && rep[i].sec-k >= 0)
					return i;
		}		
	return -1;
}
	
void extend_gapped_repeat()
{
	int i,j,k, len = rep.size();
	
	for(i =0;i<len;i++)
		if(rep[i].visited==0)
			for(j=1;j<=gap_len;j++){
				
				k = check_extend(i,j);
				if(k==-1)
					continue;
			
				//extend repeat
				rep[i].len += j+ rep[k].len-1;
				rep[i].seclen += rep[k].sec -(rep[i].sec+rep[i].seclen -1) + rep[k].seclen-1;
				rep[i].exact = 1;
				rep[k].visited = 1;
				i--;
				break;
			}
}	


void write_dna()
{
// write the dna sequence we are testing to a file. This is for debugging purposes!
    FILE *f;
	char outputfile[300];

    sprintf(outputfile, "%s.prophage.%hu.fasta", inputfile, (unsigned short) ppno);
    f = fopen(outputfile,"w");
    fprintf(f, ">pp_%d\n%s\n", ppno, dna);
    fclose(f);
}

void print_output()
{
	int i,j;
	FILE *f;
	char outputfile[300];


    sprintf(outputfile, "%s.prophage.%hu.repeatfinder", inputfile, (unsigned short) ppno);
    f = fopen(outputfile,"w");
    fprintf(f, "First repeat start\tFirst repeat end\tSecond repeat start\tSecond repeat end\tFirst len\tSecond len\tExact\tVisited\n");

	j = rep.size();

	int totalRep = 0;
	for(i=0;i<j;i++)
		if(rep[i].visited==0 && rep[i].len>= output_rep_len){
			fprintf(f,"%d\t%d\t",rep[i].fst,rep[i].fst+rep[i].len-1);
			if(rep[i].sec>-1)
				fprintf(f,"%d\t%d",rep[i].sec,rep[i].sec+rep[i].seclen-1);
			else
				fprintf(f,"%d\t%d",rep[i].sec*(-1),(rep[i].sec+rep[i].seclen-1)*(-1));
			fprintf(f, "\t%d\t%d\t%d\t%d\n", rep[i].len, rep[i].seclen, rep[i].exact, rep[i].visited);
			totalRep++;
		}
	fclose(f);
}


void run() {

	//initialize
	converter[(int) 'A']=0;
	converter[(int) 'a']=0;
	converter[(int) 'C']=1;
	converter[(int) 'c']=1;
	converter[(int) 'G']=2;
	converter[(int) 'g']=2;
	converter[(int) 'T']=3;
	converter[(int) 't']=3;
	complement[(int) 'A']='T';
	complement[(int) 'a']='t';
	complement[(int) 'C']='G';
	complement[(int) 'c']='g';
	complement[(int) 'G']='C';
	complement[(int) 'g']='c';
	complement[(int) 'T']='A';
	complement[(int) 't']='a';

	find_repeats();
	extend_repeats();

	gap_len++;
	//printf("total rep without join = %d\n",rep.size());

	if(gap_len>1)
		extend_gapped_repeat();

}


int main(int argc, char **argv)
{
	int i;

	/*if(argc<2){
		fprintf(stderr, "Command line not valid.\n -f \"fileName\" -g \"Gap Length\"\nUse 0 for no gap\n");
		exit(0);
	}*/

	for(i=1;i<argc;i++) {
		if(!strcmp(argv[i], "-f")  || !strcmp(argv[i], "-F")) {
			assert(i+1 < argc);
			strcpy(inputfile, argv[i+1]);
			i++;
		}
		else if(!strcmp(argv[i], "-g")  || !strcmp(argv[i], "-G")) {
			assert(i+1 < argc);
			sscanf(argv[i+1], "%d", &gap_len);
			i++;
		}
		/*		else if(!strcmp(argv[i], "-l")  || !strcmp(argv[i], "-L")) {
			assert(i+1 < argc);
			sscanf(argv[i+1], "%u", &INIT_DNA_LEN);
			i++;
			}*/
		else if(!strcmp(argv[i], "-l")  || !strcmp(argv[i], "-L")) {
			assert(i+1 < argc);
			sscanf(argv[i+1], "%u", &output_rep_len);
			i++;
		}
		else {
			fprintf(stderr, "Command line not valid. -f \"fileName\" \n-g \"Gap Length\" (Use 0 for not joining) \n-l \"dna length\" (Specify DNA length if its > 10 million)\n");
			exit(0);
		}
	}

	input();
    run();
    print_output();
    write_dna();
	return 0;
}

PyObject *
python_input(PyObject *self, PyObject *args) {
    /* Parse arguments */
    if(!PyArg_ParseTuple(args, "siii", &dna, &gap_len, &ppno, &debug)) {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse the arguments to python_input");
        return NULL;
    }
    dna_len = strlen(dna);
    run();

    if (debug) {
        strcpy(inputfile, "DEBUGGING_REPEATFINDER");
        print_output();
        write_dna();
    }

    // convert our vector of arrays to a Python object
    Py_ssize_t len = rep.size();
    PyObject *result = PyList_New(0);


    int totalRep = 0;
	for(int i=0;i<len;i++) {
		if(rep[i].visited==0 && rep[i].len>= output_rep_len){
            // we just use a temp variable to deal with the
            // case when second start < 0
		    int ss = rep[i].sec;
		    int se = rep[i].sec+rep[i].seclen-1;
		    if(rep[i].sec<0) {
		        ss = rep[i].sec*(-1);
		        se = (rep[i].sec+rep[i].seclen-1)*(-1);
		    }
			totalRep++;

		    PyObject *item = Py_BuildValue("{s:i, s:i,s:i,s:i,s:i}",
		    "repeat_number", totalRep,
            "first_start", rep[i].fst,
            "first_end", rep[i].fst+rep[i].len-1,
            "second_start", ss,
            "second_end", se
            );
           PyList_Append(result, item);
           }
    }

    // CRITICAL: We need to reset the repeats before we
    // call this method again from within the same run!
    for (int i=0; i<HASH_LEN; i++)
        allrepeats[i].clear();

    rep.clear();
    return result;
}


PyMODINIT_FUNC PyInit_PhiSpyRepeatFinder(void) {
    return PyModule_Create(&PhiSpyRepeatFinderModule);
}
