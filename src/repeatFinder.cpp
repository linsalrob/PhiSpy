#include <assert.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>

#include <stdlib.h>
#include <string.h>
using namespace std;

#define REPEAT_LEN                    11
#define HASH_LEN     (1<<(REPEAT_LEN*2))

unsigned INIT_DNA_LEN = 120000000; 

char inputfile[256] = "sajia.fasta";
char *dna;
int converter[128], complement[128];
int dna_len;
int output_rep_len = REPEAT_LEN;

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
int gap_len = 0;//allow gap < gap_len

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
		key = (key<<2) + converter[dna[start]];
	allrepeats[key].push_back(0);

	for(start = 1;start < dna_len-REPEAT_LEN+1; start++)
	{
		key = ((key&((1<<keylen)-1))<< 2) + converter[dna[start+REPEAT_LEN-1]];
		allrepeats[key].push_back(start);
	}

	//find reverse repeat
	key = 0;
	for(start = dna_len-1;start >dna_len-1-REPEAT_LEN; start--)
		key = (key<<2) + converter[complement[dna[start]]];
	allrepeats[key].push_back((dna_len-1)*(-1));
	
	for(start = dna_len-2;start >REPEAT_LEN-2; start--)
	{
		key= ((key&((1<<keylen)-1))<<2) + converter[complement[dna[start-REPEAT_LEN+1]]];
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
		if(dna[fst-1] == complement[dna[sec+1]])
			return;
///////////////

	k=0;
	for(i = REPEAT_LEN+fst, j = sec - REPEAT_LEN ; i<dna_len && j>-1 && i<j; i++, j--)
		if (dna[i] == complement[dna[j]])
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
		key = (key<<2) + converter[dna[i]];
	for(j=0;j<allrepeats[key].size();j++)
		if(allrepeats[key][j]<0)
			find_maxlen_rev(0,allrepeats[key][j]);
		else
			find_maxlen(0,allrepeats[key][j]);

	for(i =1;i<dna_len-REPEAT_LEN+1;i++){
		key = ((key&((1<<keylen)-1))<< 2) + converter[dna[i+REPEAT_LEN-1]];
		for(j=0;j<allrepeats[key].size();j++)
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
			tail = mid -1;
	}
	i = mid-1;
	while(rep[mid].fst == rep[i].fst && i >-1)
		i--;
	i++;
////
	k = rep[fst].sec+rep[fst].seclen -1;
	
	for(;rep[i].fst == j;i++)
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

void print_output()
{
	int i,j,k;
	FILE *f;
	char outputfile[256];

        strcpy (outputfile,inputfile);
        strcat (outputfile,".repeatfinder");
        f = fopen(outputfile,"w");

	j = rep.size();

	int totalRep = 0;
	for(i=0;i<j;i++)
		if(rep[i].visited==0 && rep[i].len>= output_rep_len){
		        
		  ////print
			fprintf(f,"%d\t%d\t",rep[i].fst,rep[i].fst+rep[i].len-1);
			
			if(rep[i].sec>-1){
				fprintf(f,"%d\t%d\n",rep[i].sec,rep[i].sec+rep[i].seclen-1);	
				/*
				for(k = rep[i].fst-1; k<rep[i].fst+rep[i].len-1;k++) fprintf(f,"%c",dna[k]);
				fprintf(f,"\t");
			
				for(k = rep[i].sec-1; k<rep[i].sec+rep[i].seclen-1;k++) 
					fprintf(f,"%c",dna[k]);
					fprintf(f,"\tsame direction\t");*/
			}
			else{
				fprintf(f,"%d\t%d\n",rep[i].sec*(-1),(rep[i].sec+rep[i].seclen-1)*(-1));
			        /* 
				for(k = rep[i].fst-1; k<rep[i].fst+rep[i].len-1;k++) fprintf(f,"%c",dna[k]);
				fprintf(f,"\t");
			
				for(k = rep[i].sec+rep[i].seclen; k>rep[i].sec;k--) 
					fprintf(f,"%c",dna[k*(-1)]);
				fprintf(f,"\treverse complement\t");*/
			}
			/*if(rep[i].exact==0)
				fprintf(f,"exact\n");
			else
				fprintf(f,"inexact\n");
			*/			
			totalRep++;
		}
	//printf("total repeat (join <= %d nt) = %d\n",gap_len-1,totalRep);
	fclose(f);
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

	//initialize
	converter['A']=0;
	converter['a']=0;
	converter['C']=1;
	converter['c']=1;
	converter['G']=2;
	converter['g']=2;
	converter['T']=3;
	converter['t']=3;
	complement['A']='T';
	complement['a']='t';
	complement['C']='G';
	complement['c']='g';
	complement['G']='C';
	complement['g']='c';
	complement['T']='A';
	complement['t']='a';
	
	input();
	find_repeats();
	extend_repeats();

	gap_len++;
	//printf("total rep without join = %d\n",rep.size());

	if(gap_len>1)
		extend_gapped_repeat();
		
	print_output();
	return 0;
}
