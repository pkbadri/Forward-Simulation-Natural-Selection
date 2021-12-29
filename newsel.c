//Debugged on 08/23/2011
#include "mtrand.h"
#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include <fstream.h>
#include <time.h>
#include <assert.h> 
#define       ZERO   0          //Constant. Do not change.

//FIX THE NUMBER OF SITES UNDER SELECTION HERE.
static double N;                //Number of sites under selection. Choose a small value (e.g. 1-200).

//SITES ARE ASSUMED TO BE BIALLELIC AND THERE ARE 3 POSSIBLE GENOTYPES AT A SITE 00, 10 and 11. 
//SET THESE ARRAYS in newselinput.txt. FITNESS FOR EACH GENOTYPE IS 1 + SELECTION COEFFICIENT.
static int  position[200];      //Positions of selected sites in basepairs in ascending order. [0 - SEQLENGTH-1]
static double      S[200];      //Strength of selection per site for the selected allele homozygote.(11)
static double      H[200];      //Strength of selection per site for heterozygotes.(01 or 10)
static double     NU[200];      //Strength of selection per site for the alternate allele homozygote.(00)

//DO NOT CHANGE VARIABLES BELOW THIS.
static int         SSIZE;       //Number of samples to output from the final population. Choose <= POPSIZE/2 for chromosomes from different individuals.
static int       POPSIZE;       //Total number of chromsomes in the diploid population.
static int     SEQLENGTH;       //Length of the DNA sequence in basepairs. Choose >= DELETE*mu*POPSIZE. 
static double         re;       //Per-generation per-sequence rate of recombination.
static double         mu;       //Per-generation per-sequence rate of mutation.    
static double          s;       //Probability of self-fertilization.
static int           GEN;       //Number of generations.
static int        DELETE;       //Intervals after which fixed mutations are removed.  
static int       FITNESS;       //Fitness effects type


unsigned long init[4] = {0x123, 0x234, 0x345, 0x456},length = 4;
unsigned long duplicate[4];
MTRand drand;//double in [0, 1) generator, already init  
MTRand mt;

static inline double  minimum(double i,double j){return(i > j) ? j:i;}
static inline double  maximum(double i,double j){return(i > j) ? i:j;}
static inline int poisson(float mean){int count=0;double q,bound;bound = exp(-mean);
for(q=1;q>=bound;q*=drand()){++count;}return(count-1);}
static int *indicator,*multhit,*selected,*tmpp;

//Ancestry information for a generation. 
class next{
public:
double *np;int *n, *brk;
};

//Arrays representing chromosomes
class member{
public:
int *q;//Choose a suitable value depending on the expected number of mutations per individual
};static member *temp;

//Pointers to chromosome arrays
struct node{
member *chr;
};static node *odd,*even; 


int i,j,k,l,m,n,o,r,trp,trk,start,end,ind,mut,RUNS;
static int *n2count,*n2r,*n3count,*n3r,*n4count,*n4r,*n5count,*n5r;
static int *n6count,*n6r,*n7count,*n7r,*n8count,*n8r,*n9count,*n9r;
static int *current, *size, *siz, *mcount,*recount;
static int *temp1, *temp2, *indic, flag1,flag2;

double V,p,rec,x,y,z;
short  (*selodd)[200], (*seleven)[200], (*fsample)[200];
double *efit, *ofit, maxfit;
next   *gen1,*gen2,*gen3,*gen4,*gen5,*gen6,*gen7,*gen8,*gen9,*gtemp;

 
main(int argc, char *argv[]){
SSIZE = atoi(argv[2]);POPSIZE = atoi(argv[4]);SEQLENGTH = atoi(argv[6]);re = atof(argv[8]);
mu = atof(argv[10]);s = atof(argv[12]);GEN = atoi(argv[14]);DELETE = atoi(argv[16]);FITNESS = atoi(argv[20]);


if(POPSIZE % 2 == 1){cout<<"Error. pop must be even number.\n"; return(0);}

 
n2count = new int[POPSIZE];n2r = new int[POPSIZE];
n3count = new int[POPSIZE];n3r = new int[POPSIZE];
n4count = new int[POPSIZE];n4r = new int[POPSIZE];
n5count = new int[POPSIZE];n5r = new int[POPSIZE];
n6count = new int[POPSIZE];n6r = new int[POPSIZE];
n7count = new int[POPSIZE];n7r = new int[POPSIZE];
n8count = new int[POPSIZE];n8r = new int[POPSIZE];
n9count = new int[POPSIZE];n9r = new int[POPSIZE];
current = new int[POPSIZE];size = new int[POPSIZE];siz = new int[POPSIZE];
temp1 = new int[POPSIZE];temp2 = new int[POPSIZE];
recount = new int[POPSIZE];indic = new int [POPSIZE];
mcount = new int[POPSIZE];

selodd   = new short[POPSIZE][200]; seleven = new short[POPSIZE][200];
fsample  = new short[POPSIZE][200]; 
efit    = new double[POPSIZE]; ofit = new double[POPSIZE];


ifstream input("newselinput.txt",ios::in);
while(!input.eof()){
input>>N;if(N > 200){cout<<"First line in newselinput.txt should be <= 200.\n";return(0);}
for(i=0;i<N;i++){input>>position[i];}
for(i=0;i<N;i++){input>>H[i];  if(FITNESS == 2){H[i]  +=  1.0;}}
for(i=0;i<N;i++){input>>S[i];  if(FITNESS == 2){S[i]  +=  1.0;}}
for(i=0;i<N;i++){input>>NU[i]; if(FITNESS == 2){NU[i] +=  1.0;}}
break;}
input.close();

selected = new int[SEQLENGTH];memset(&selected[0],'\0',4*SEQLENGTH);

//DO NOT CHANGE LINES BELOW THIS.
for(i=0;i<N;i++){selected[position[i]] = 1;}

rec = 1 - exp(-re);

gen1 = new next; gen1->np = new double[POPSIZE]; gen1->n = new int[POPSIZE]; gen1->brk = new int[POPSIZE];
gen2 = new next; gen2->np = new double[POPSIZE]; gen2->n = new int[POPSIZE]; gen2->brk = new int[POPSIZE];
gen3 = new next; gen3->np = new double[POPSIZE]; gen3->n = new int[POPSIZE]; gen3->brk = new int[POPSIZE];
gen4 = new next; gen4->np = new double[POPSIZE]; gen4->n = new int[POPSIZE]; gen4->brk = new int[POPSIZE];
gen5 = new next; gen5->np = new double[POPSIZE]; gen5->n = new int[POPSIZE]; gen5->brk = new int[POPSIZE];
gen6 = new next; gen6->np = new double[POPSIZE]; gen6->n = new int[POPSIZE]; gen6->brk = new int[POPSIZE];
gen7 = new next; gen7->np = new double[POPSIZE]; gen7->n = new int[POPSIZE]; gen7->brk = new int[POPSIZE];
gen8 = new next; gen8->np = new double[POPSIZE]; gen8->n = new int[POPSIZE]; gen8->brk = new int[POPSIZE];
gen9 = new next; gen9->np = new double[POPSIZE]; gen9->n = new int[POPSIZE]; gen9->brk = new int[POPSIZE];

gtemp = new next;gtemp->np = new double[POPSIZE];gtemp->n = new int[POPSIZE];gtemp->brk = new int[POPSIZE];
temp = new member;odd = new node[POPSIZE];even = new node[POPSIZE];
indicator = new int[SEQLENGTH]; multhit = new int[SEQLENGTH];
V = (N*mu/SEQLENGTH); mu -= V;


for(RUNS=1; RUNS<=atoi(argv[18]); RUNS++){
 
for(i=0;i<length;i++){init[i] = time(NULL); duplicate[i] = init[i];} mt.seed(init,length);

int array_size, finished, flag;  array_size = 1500 + int(2.5*POPSIZE*mu);


finished = 0;
while (finished == 0){

array_size += 100;
memset(&size[0],   '\0', 4*POPSIZE);
memset(&siz[0],    '\0', 4*POPSIZE);
memset(&multhit[0],'\0', 4*SEQLENGTH);

flag = 0;
for(i=0;i<POPSIZE;i++){odd[i].chr = new member; odd[i].chr->q = new int[array_size]; even[i].chr = new member;even[i].chr->q = new int[array_size];}
for(i=0;i<length;i++){init[i] = duplicate[i];} mt.seed(init, length);

//Initialize the size of chromosome arrays and arrays
for(i=0;i<POPSIZE;i++){
for(j=0; j<N; j++){seleven[i][j] = 0;} 
efit[i] = 1;}


//Create generation 1
m = 0; maxfit = 0; 
while(m < POPSIZE){

//Choose chromosomes according to fitness
gen1->np[m] = drand(); gen1->np[m+1] = drand(); 
if(gen1->np[m] < rec){gen1->brk[m] = int((SEQLENGTH - 1)*drand());}  if(gen1->np[m+1] < rec){gen1->brk[m+1]  =  int((SEQLENGTH  -  1)*drand());}

x = drand();
//If not created by selfing
if(x >= s){
gen1->n[m]   = int(POPSIZE*drand());  z = drand(); while(z >= efit[gen1->n[m]])  {gen1->n[m]   = int(POPSIZE*drand());  z = drand();} 
gen1->n[m+1] = int(POPSIZE*drand());  z = drand(); while(z >= efit[gen1->n[m+1]]){gen1->n[m+1] = int(POPSIZE*drand());  z = drand();}
}

//If individual is created by selfing
else{
y = drand();
if(y <= 0.50){
gen1->n[m] = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen1->n[m]]){gen1->n[m] = int(POPSIZE*drand());z=drand();}
gen1->n[m+1] = gen1->n[m];}

else{
gen1->n[m]  =  int(POPSIZE*drand());  z  =  drand();while(z >= efit[gen1->n[m]]){gen1->n[m]=int(POPSIZE*drand());z = drand();}
if((gen1->n[m])%2  == 0){gen1->n[m+1] = (gen1->n[m]) + 1;}  if((gen1->n[m])%2 == 1){gen1->n[m+1]  =  (gen1->n[m])  -  1;}
}}

 
//Create the next generation
j = gen1->n[m];

//If the newly chosen chromosome is a recombinant
if(gen1->np[m] < rec){
if((j&1) == 0){k = j+1;}  else{k = j-1;}
l = 0;while(position[l] <= gen1->brk[m]){selodd[m][l]=seleven[j][l];++l;if(l==N){break;}}
memcpy(&selodd[m][l],&seleven[k][l],2*(N - l));}

//If newly chosen chromosome is a non-recombinant
else{memcpy(&selodd[m][0],  &seleven[j][0], 2*N);}  

//Create the next generation
j = gen1->n[m+1];

//If the newly chosen chromosome is a recombinant
if(gen1->np[m+1] < rec){
if((j&1)==0){k = j+1;}else{k = j-1;}
l = 0; while(position[l] <= gen1->brk[m+1]){selodd[m+1][l] = seleven[j][l]; ++l; if(l == N){break;}}
memcpy(&selodd[m+1][l], &seleven[k][l], 2*(N - l));
}

//If newly chosen chromosome is a non-recombinant
else{memcpy(&selodd[m+1][0], &seleven[j][0],2*N);}

//Mutate selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m][j]   == 0){selodd[m][j]   = 1;continue;}if(selodd[m][j]   == 1){selodd[m][j]   = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m+1][j] == 0){selodd[m+1][j] = 1;continue;}if(selodd[m+1][j] == 1){selodd[m+1][j] = 0;continue;}}

//Calculate fitness for each individual in the population
ofit[m] = 1; ofit[m+1] = 1;
for(l=0;l<N;l++){
if(selodd[m][l]      != selodd[m+1][l]     ){if(FITNESS == 1){ofit[m] +=  H[l]; ofit[m+1] +=  H[l];} else{ofit[m] *=  H[l];  ofit[m+1] *= H[l]; }}
if(selodd[m][l] == 1 && selodd[m+1][l] == 1){if(FITNESS == 1){ofit[m] +=  S[l]; ofit[m+1] +=  S[l];} else{ofit[m] *=  S[l];  ofit[m+1] *= S[l]; }}
if(selodd[m][l] == 0 && selodd[m+1][l] == 0){if(FITNESS == 1){ofit[m] += NU[l]; ofit[m+1] += NU[l];} else{ofit[m] *=  NU[l]; ofit[m+1] *= NU[l];}}
}

if(ofit[m] > maxfit){maxfit = ofit[m];}
m = m + 2;}

for(m = 0; m < POPSIZE; m++){ofit[m] /= maxfit;}

 


//Create generation 2
m = 0; maxfit = 0;
while(m  <  POPSIZE){
//Choose chromosomes according to fitness
gen2->np[m] = drand(); gen2->np[m+1] = drand();
if(gen2->np[m] < rec){gen2->brk[m] = int((SEQLENGTH-1)*drand());}if(gen2->np[m+1] < rec){gen2->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
//If not created by selfing
if(x >= s){
gen2->n[m]   = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen2->n[m]]  ){gen2->n[m]   = int(POPSIZE*drand()); z = drand();}
gen2->n[m+1] = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen2->n[m+1]]){gen2->n[m+1] = int(POPSIZE*drand()); z = drand();}
}

//If created by selfing
else{
y = drand();
if(y <= 0.50){
gen2->n[m] = int(POPSIZE*drand());z = drand(); while(z >= ofit[gen2->n[m]]){gen2->n[m]  =  int(POPSIZE*drand()); z = drand();}
gen2->n[m+1]= gen2->n[m];
}

else{
gen2->n[m] = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen2->n[m]]){gen2->n[m] = int(POPSIZE*drand()); z = drand();}
if((gen2->n[m])%2 == 0){gen2->n[m+1]  =  (gen2->n[m])+1;}  if((gen2->n[m])%2  ==  1){gen2->n[m+1]  =  (gen2->n[m])-1;}}
}

//Create the next generation
j = gen2->n[m];if(gen2->np[m] < rec){if((j&1)==0){k = j + 1;}else{k = j - 1;}
l=0;while(position[l] <= gen2->brk[m]){seleven[m][l] = selodd[j][l];++l;if(l==N){break;}}memcpy(&seleven[m][l],&selodd[k][l],2*(N-l));}
else{memcpy(&seleven[m][0],&selodd[j][0],2*N);}
j = gen2->n[m+1];if(gen2->np[m+1] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen2->brk[m+1]){seleven[m+1][l]=selodd[j][l];++l;if(l==N){break;}}memcpy(&seleven[m+1][l],&selodd[k][l],2*(N-l));}
else{memcpy(&seleven[m+1][0],&selodd[j][0],2*N);}


//Mutate the selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(seleven[m][j]==0){seleven[m][j] = 1;continue;}if(seleven[m][j]==1){seleven[m][j] = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(seleven[m+1][j]==0){seleven[m+1][j] = 1;continue;}if(seleven[m+1][j]==1){seleven[m+1][j] = 0;continue;}}

//Calculate the fitness of each individual
efit[m] = 1; efit[m+1] = 1;
for(l=0;l<N;l++){
if(seleven[m][l]      !=     seleven[m+1][l]){if(FITNESS == 1){efit[m] += H[l]; efit[m+1] += H[l];}  else{efit[m] *= H[l];  efit[m+1]  *= H[l];}}
if(seleven[m][l] == 1 && seleven[m+1][l] ==1){if(FITNESS == 1){efit[m] += S[l]; efit[m+1] += S[l];}  else{efit[m] *= S[l];  efit[m+1]  *= S[l];}}
if(seleven[m][l] == 0 && seleven[m+1][l] ==0){if(FITNESS == 1){efit[m] += NU[l];efit[m+1] += NU[l];} else{efit[m] *= NU[l]; efit[m+1]  *= NU[l];}}
}

if(efit[m] > maxfit){maxfit = efit[m];}
m = m + 2;}

for(m=0;m<POPSIZE;m++){efit[m] /= maxfit;}




//Create generation 3 
m = 0; maxfit = 0;
while(m < POPSIZE){

//Choose chromosomes
gen3->np[m] = drand();gen3->np[m+1] = drand();x=drand();
if(gen3->np[m] < rec){gen3->brk[m] = int((SEQLENGTH-1)*drand());}if(gen3->np[m+1] < rec){gen3->brk[m+1] = int((SEQLENGTH-1)*drand());}
if(x>=s){gen3->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen3->n[m]]){gen3->n[m]=int(POPSIZE*drand());z=drand();}
gen3->n[m+1]=int(POPSIZE*drand());z=drand();while(z >= efit[gen3->n[m+1]]){gen3->n[m+1]=int(POPSIZE*drand());z=drand();}}
else{y=drand();if(y<=0.50){gen3->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen3->n[m]]){gen3->n[m]=int(POPSIZE*drand());z=drand();}
gen3->n[m+1]= gen3->n[m];}else{gen3->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen3->n[m]]){gen3->n[m]=int(POPSIZE*drand());z=drand();}
if((gen3->n[m])%2==0){gen3->n[m+1]=(gen3->n[m])+1;}if((gen3->n[m])%2==1){gen3->n[m+1]=(gen3->n[m])-1;}}}

//Create next generation
j = gen3->n[m];if(gen3->np[m] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen3->brk[m]){selodd[m][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m][0],&seleven[j][0],2*N);}
j = gen3->n[m+1];if(gen3->np[m+1] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen3->brk[m+1]){selodd[m+1][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m+1][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m+1][0],&seleven[j][0],2*N);}


//Mutate selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m][j]==0){selodd[m][j] = 1;continue;}if(selodd[m][j]==1){selodd[m][j] = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m+1][j]==0){selodd[m+1][j] = 1;continue;}if(selodd[m+1][j]==1){selodd[m+1][j] = 0;continue;}}


//Calculate fitness
ofit[m] = 1;ofit[m+1] = 1;
for(l=0;l<N;l++){
if(selodd[m][l]    != selodd[m+1][l]   ){if(FITNESS == 1){ofit[m]+= H[l];ofit[m+1]+= H[l];}else{ofit[m] *= H[l];ofit[m+1]*= H[l];}}
if(selodd[m][l]==1 && selodd[m+1][l]==1){if(FITNESS == 1){ofit[m]+= S[l];ofit[m+1]+= S[l];}else{ofit[m] *= S[l];ofit[m+1]*= S[l];}}
if(selodd[m][l]==0 && selodd[m+1][l]==0){if(FITNESS == 1){ofit[m]+= NU[l];ofit[m+1]+= NU[l];}else{ofit[m] *= NU[l];ofit[m+1]*= NU[l];}}
}

if(ofit[m] > maxfit){maxfit = ofit[m];}
m = m + 2;}

for(m=0;m<POPSIZE;m++){ofit[m]  /=  maxfit;}




//Create generation 4
m = 0; maxfit = 0;
while(m < POPSIZE){

//Choose chromosomes
gen4->np[m] = drand(); gen4->np[m+1] = drand();
if(gen4->np[m] < rec){gen4->brk[m] = int((SEQLENGTH-1)*drand());}if(gen4->np[m+1] < rec){gen4->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
if(x >= s){
gen4->n[m]   = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen4->n[m]]  ){gen4->n[m]   = int(POPSIZE*drand()); z = drand();}
gen4->n[m+1] = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen4->n[m+1]]){gen4->n[m+1] = int(POPSIZE*drand()); z = drand();}
}

else{
y = drand();
if(y <= 0.50){
gen4->n[m] = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen4->n[m]]){gen4->n[m]  =  int(POPSIZE*drand());  z = drand();}
gen4->n[m+1]= gen4->n[m];
}

else{
gen4->n[m] = int(POPSIZE*drand());z = drand();while(z >= ofit[gen4->n[m]]){gen4->n[m]=int(POPSIZE*drand()); z = drand();}
if((gen4->n[m])%2 == 0){gen4->n[m+1] = (gen4->n[m]) + 1;} if((gen4->n[m])%2  ==  1){gen4->n[m+1] = (gen4->n[m])-1;}
}
}


//Create next generation
j = gen4->n[m];   if(gen4->np[m]   < rec){if((j&1) == 0){k = j + 1;} else{k = j - 1;}
l = 0; while(position[l] <= gen4->brk[m]){seleven[m][l]  =  selodd[j][l]; ++l; if(l == N){break;}} memcpy(&seleven[m][l], &selodd[k][l], 2*(N-l));}
else{memcpy(&seleven[m][0], &selodd[j][0], 2*N);}

j = gen4->n[m+1]; 
if(gen4->np[m+1] < rec){
if((j&1) == 0){k = j + 1;} else{k = j - 1;}
 l = 0; while(position[l] <= gen4->brk[m+1]){seleven[m+1][l] = selodd[j][l]; ++l; if(l == N){break;}} memcpy(&seleven[m+1][l],&selodd[k][l],2*(N-l));
}
else{memcpy(&seleven[m+1][0],  &selodd[j][0],  2*N);}

//Mutate the selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(seleven[m][j]   == 0){seleven[m][j]   = 1;continue;}if(seleven[m][j]   == 1){seleven[m][j]   = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(seleven[m+1][j] == 0){seleven[m+1][j] = 1;continue;}if(seleven[m+1][j] == 1){seleven[m+1][j] = 0;continue;}}

//Calculate fitness
efit[m] = 1; efit[m+1] = 1;
for(l=0;l<N;l++){
if(seleven[m][l]       !=  seleven[m+1][l]     ){if(FITNESS == 1){efit[m] +=  H[l]; efit[m+1] += H[l];  }else{efit[m] *= H[l];  efit[m+1] *= H[l];  }}
if(seleven[m][l] == 1  &&  seleven[m+1][l] == 1){if(FITNESS == 1){efit[m] +=  S[l]; efit[m+1] += S[l];  }else{efit[m] *= S[l];  efit[m+1] *= S[l];  }}
if(seleven[m][l] == 0  &&  seleven[m+1][l] == 0){if(FITNESS == 1){efit[m] += NU[l]; efit[m+1] += NU[l]; }else{efit[m] *= NU[l]; efit[m+1] *= NU[l]; }}
}

if(efit[m] > maxfit){maxfit = efit[m];}
m = m + 2;}

for(m = 0; m < POPSIZE; m++){efit[m] /= maxfit;}




//Create  generation 5
m = 0; maxfit = 0;
while(m < POPSIZE){

//Choose chromosomes
gen5->np[m] = drand();gen5->np[m+1] = drand();x=drand();
if(gen5->np[m] < rec){gen5->brk[m] = int((SEQLENGTH-1)*drand());}if(gen5->np[m+1] < rec){gen5->brk[m+1] = int((SEQLENGTH-1)*drand());}

if(x>=s){gen5->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen5->n[m]]){gen5->n[m]=int(POPSIZE*drand());z=drand();}
gen5->n[m+1]=int(POPSIZE*drand());z=drand();while(z >= efit[gen5->n[m+1]]){gen5->n[m+1]=int(POPSIZE*drand());z=drand();}}
else{y=drand();if(y<=0.50){gen5->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen5->n[m]]){gen5->n[m]=int(POPSIZE*drand());z=drand();}
gen5->n[m+1]= gen5->n[m];}else{gen5->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen5->n[m]]){gen5->n[m]=int(POPSIZE*drand());z=drand();}
if((gen5->n[m])%2==0){gen5->n[m+1]=(gen5->n[m])+1;}if((gen5->n[m])%2==1){gen5->n[m+1]=(gen5->n[m])-1;}}}


//Create next generation
j = gen5->n[m];if(gen5->np[m] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen5->brk[m]){selodd[m][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m][0],&seleven[j][0],2*N);}
j = gen5->n[m+1];if(gen5->np[m+1] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen5->brk[m+1]){selodd[m+1][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m+1][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m+1][0],&seleven[j][0],2*N);}


//Mutate selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m][j]==0){selodd[m][j] = 1;continue;}if(selodd[m][j]==1){selodd[m][j] = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m+1][j]==0){selodd[m+1][j] = 1;continue;}if(selodd[m+1][j]==1){selodd[m+1][j] = 0;continue;}}

//Calculate fitness
ofit[m]=1;ofit[m+1]=1;
for(l=0;l<N;l++){
if(selodd[m][l]!=selodd[m+1][l]){if(FITNESS==1){ofit[m]+= H[l];ofit[m+1]+= H[l];}else{ofit[m] *= H[l];ofit[m+1]*= H[l];}}
if(selodd[m][l]==1 && selodd[m+1][l]==1){if(FITNESS==1){ofit[m]+= S[l];ofit[m+1]+= S[l];}else{ofit[m] *= S[l];ofit[m+1]*= S[l];}}
if(selodd[m][l]==0 && selodd[m+1][l]==0){if(FITNESS==1){ofit[m]+= NU[l];ofit[m+1]+= NU[l];}else{ofit[m] *= NU[l];ofit[m+1]*= NU[l];}}
}

if(ofit[m] > maxfit){maxfit = ofit[m];}
m = m + 2;}

for(m=0; m < POPSIZE; m++){ofit[m] /= maxfit;}




//Create generation 6 
m = 0; maxfit = 0;
while(m < POPSIZE){

//Choose chromosomes
gen6->np[m] = drand();gen6->np[m+1] = drand();
if(gen6->np[m] < rec){gen6->brk[m] = int((SEQLENGTH-1)*drand());}if(gen6->np[m+1] < rec){gen6->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
if(x >= s){
gen6->n[m]   = int(POPSIZE*drand()); z = drand();while(z >= ofit[gen6->n[m]]  ){gen6->n[m]   = int(POPSIZE*drand()); z = drand();}
gen6->n[m+1] = int(POPSIZE*drand()); z = drand();while(z >= ofit[gen6->n[m+1]]){gen6->n[m+1] = int(POPSIZE*drand()); z = drand();}
}

else{
y = drand();if(y <= 0.50){gen6->n[m] = int(POPSIZE*drand());   z = drand();while(z >= ofit[gen6->n[m]]){gen6->n[m]=int(POPSIZE*drand());z=drand();}
gen6->n[m+1]= gen6->n[m];}else{gen6->n[m]=int(POPSIZE*drand());z = drand();while(z >= ofit[gen6->n[m]]){gen6->n[m]=int(POPSIZE*drand());z=drand();}
if((gen6->n[m])%2 == 0){gen6->n[m+1] = (gen6->n[m])+1;}if((gen6->n[m])%2 == 1){gen6->n[m+1] = (gen6->n[m])-1;}}
}


//Create next generation
j = gen6->n[m];
if(gen6->np[m] < rec){
if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0;while(position[l] <= gen6->brk[m]  ){seleven[m][l]  = selodd[j][l];   ++l; if(l == N){break;}} memcpy(&seleven[m][l],&selodd[k][l],2*(N-l));
}
else{memcpy(&seleven[m][0],  &selodd[j][0],  2*N);}

j = gen6->n[m+1];
if(gen6->np[m+1] < rec){
if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0;while(position[l] <= gen6->brk[m+1]){seleven[m+1][l] =  selodd[j][l]; ++l; if(l == N){break;}} memcpy(&seleven[m+1][l],&selodd[k][l],2*(N-l));
}
else{memcpy(&seleven[m+1][0], &selodd[j][0], 2*N);}

//Mutate at selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(seleven[m][j]   == 0){seleven[m][j]   = 1;continue;}if(seleven[m][j]   == 1){seleven[m][j]   = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(seleven[m+1][j] == 0){seleven[m+1][j] = 1;continue;}if(seleven[m+1][j] == 1){seleven[m+1][j] = 0;continue;}}

//Calculate fitness
efit[m] = 1; efit[m+1] = 1;
for(l=0;l<N; l++){
if(seleven[m][l]       !=  seleven[m+1][l]     ){if(FITNESS == 1){efit[m] += H[l]; efit[m+1] += H[l];}  else{efit[m] *= H[l]; efit[m+1] *= H[l];}}
if(seleven[m][l] == 1  &&  seleven[m+1][l] == 1){if(FITNESS == 1){efit[m] += S[l]; efit[m+1] += S[l];}  else{efit[m] *= S[l]; efit[m+1] *= S[l];}}
if(seleven[m][l] == 0  &&  seleven[m+1][l] == 0){if(FITNESS == 1){efit[m] += NU[l];efit[m+1] += NU[l];} else{efit[m] *= NU[l];efit[m+1] *= NU[l];}}
}

if(efit[m] > maxfit){maxfit = efit[m];}
m = m + 2;}

for(m=0; m < POPSIZE; m++){efit[m] /= maxfit;}




//Create generation 7
m = 0; maxfit = 0;
while(m < POPSIZE){

//Choose chromosomes
gen7->np[m] = drand();gen7->np[m+1] = drand();
if(gen7->np[m] < rec){gen7->brk[m] = int((SEQLENGTH-1)*drand());}if(gen7->np[m+1] < rec){gen7->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
if(x >= s){
gen7->n[m]    = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen7->n[m]]  ){gen7->n[m]   = int(POPSIZE*drand()); z = drand();}
gen7->n[m+1]  = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen7->n[m+1]]){gen7->n[m+1] = int(POPSIZE*drand()); z = drand();}
}
else{
y = drand();
if(y <= 0.50){gen7->n[m] = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen7->n[m]]){gen7->n[m] = int(POPSIZE*drand()); z = drand();}
gen7->n[m+1] = gen7->n[m];} else{gen7->n[m] = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen7->n[m]]){gen7->n[m] = int(POPSIZE*drand()); z = drand();}
if((gen7->n[m])%2 == 0){gen7->n[m+1] = (gen7->n[m]) + 1;}if((gen7->n[m])%2 == 1){gen7->n[m+1] = (gen7->n[m])-1;}}
}


//Create next generation
j = gen7->n[m];if(gen7->np[m] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen7->brk[m]){selodd[m][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m][0],&seleven[j][0],2*N);}
j = gen7->n[m+1];if(gen7->np[m+1] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen7->brk[m+1]){selodd[m+1][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m+1][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m+1][0],&seleven[j][0],2*N);}

//Mutate at selected sites
k = poisson(V); for(l=0;l<k;l++){j = int(N*drand()); if(selodd[m][j]   == 0){selodd[m][j]   = 1;continue;}if(selodd[m][j]   == 1){selodd[m][j]   = 0;continue;}}
k = poisson(V); for(l=0;l<k;l++){j = int(N*drand()); if(selodd[m+1][j] == 0){selodd[m+1][j] = 1;continue;}if(selodd[m+1][j] == 1){selodd[m+1][j] = 0;continue;}}

//Calculate fitness
ofit[m] = 1; ofit[m+1] = 1;
for(l=0;l<N;l++){
if(selodd[m][l]    !=    selodd[m+1][l]){if(FITNESS==1){ofit[m]+= H[l];ofit[m+1]+= H[l];}else{ofit[m] *= H[l];ofit[m+1]*= H[l];}}
if(selodd[m][l]==1 && selodd[m+1][l]==1){if(FITNESS==1){ofit[m]+= S[l];ofit[m+1]+= S[l];}else{ofit[m] *= S[l];ofit[m+1]*= S[l];}}
if(selodd[m][l]==0 && selodd[m+1][l]==0){if(FITNESS==1){ofit[m]+= NU[l];ofit[m+1]+= NU[l];}else{ofit[m] *= NU[l];ofit[m+1]*= NU[l];}}
}

if(ofit[m] > maxfit){maxfit = ofit[m];}
m = m + 2;}

for(m=0;m<POPSIZE;m++){ofit[m] /= maxfit;}




//Create generation 8
m = 0; maxfit = 0;
while(m < POPSIZE){

//Choose chromosomes
gen8->np[m] = drand(); gen8->np[m+1] = drand();
if(gen8->np[m] < rec){gen8->brk[m] = int((SEQLENGTH-1)*drand());}if(gen8->np[m+1] < rec){gen8->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
if(x >= s){
gen8->n[m]   = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen8->n[m]]  ){gen8->n[m]   = int(POPSIZE*drand()); z = drand();}
gen8->n[m+1] = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen8->n[m+1]]){gen8->n[m+1] = int(POPSIZE*drand()); z = drand();}
}
else{
y = drand();
if(y <= 0.50){gen8->n[m]=int(POPSIZE*drand());z=drand();while(z >= ofit[gen8->n[m]]){gen8->n[m]=int(POPSIZE*drand()); z = drand();}
gen8->n[m+1]= gen8->n[m];}else{gen8->n[m]=int(POPSIZE*drand());z=drand();while(z >= ofit[gen8->n[m]]){gen8->n[m]=int(POPSIZE*drand()); z = drand();}
if((gen8->n[m])%2 == 0){gen8->n[m+1]=(gen8->n[m])+1;}if((gen8->n[m])%2==1){gen8->n[m+1]=(gen8->n[m])-1;}}
}

//Create next generation
j = gen8->n[m];if(gen8->np[m] < rec){if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0;while(position[l] <= gen8->brk[m]){seleven[m][l]=selodd[j][l];++l;if(l==N){break;}}memcpy(&seleven[m][l],&selodd[k][l],2*(N-l));}
else{memcpy(&seleven[m][0],&selodd[j][0],2*N);}
j = gen8->n[m+1];if(gen8->np[m+1] < rec){if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0;while(position[l] <= gen8->brk[m+1]){seleven[m+1][l]=selodd[j][l];++l;if(l==N){break;}}memcpy(&seleven[m+1][l],&selodd[k][l],2*(N-l));}
else{memcpy(&seleven[m+1][0],&selodd[j][0],2*N);}

//Mutate the selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand()); if(seleven[m][j]   == 0){seleven[m][j]   = 1;continue;}if(seleven[m][j]   == 1){seleven[m][j]   = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand()); if(seleven[m+1][j] == 0){seleven[m+1][j] = 1;continue;}if(seleven[m+1][j] == 1){seleven[m+1][j] = 0;continue;}}

//Calculate fitness
efit[m] = 1; efit[m+1] = 1;
for(l=0;l<N;l++){
if(seleven[m][l] !=    seleven[m+1][l]){if(FITNESS==1){efit[m]+= H[l];efit[m+1]+= H[l];}else{efit[m] *= H[l];efit[m+1]*= H[l];}}
if(seleven[m][l]==1 && seleven[m+1][l]==1){if(FITNESS==1){efit[m]+= S[l];efit[m+1]+= S[l];}else{efit[m] *= S[l];efit[m+1]*= S[l];}}
if(seleven[m][l]==0 && seleven[m+1][l]==0){if(FITNESS==1){efit[m]+= NU[l];efit[m+1]+= NU[l];}else{efit[m] *= NU[l];efit[m+1]*= NU[l];}}
}

if(efit[m] > maxfit){maxfit = efit[m];}
m = m + 2;}

for(m=0;m<POPSIZE;m++){efit[m] /= maxfit;}




//Create generation 9
m = 0; maxfit = 0;
while(m < POPSIZE){
//Choose chromosomes
gen9->np[m] = drand(); gen9->np[m+1] = drand();
if(gen9->np[m] < rec){gen9->brk[m] = int((SEQLENGTH-1)*drand());}if(gen9->np[m+1] < rec){gen9->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
if(x >= s){gen9->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen9->n[m]]){gen9->n[m]=int(POPSIZE*drand());z=drand();}
gen9->n[m+1]=int(POPSIZE*drand());z=drand();while(z >= efit[gen9->n[m+1]]){gen9->n[m+1]=int(POPSIZE*drand());z=drand();}}
else{
y = drand();
if(y <= 0.50){
gen9->n[m]   = int(POPSIZE*drand());z = drand();while(z >= efit[gen9->n[m]]){gen9->n[m] = int(POPSIZE*drand()); z = drand();}
gen9->n[m+1] = gen9->n[m];
}
else{
gen9->n[m] = int(POPSIZE*drand()); z = drand();while(z >= efit[gen9->n[m]]){gen9->n[m] = int(POPSIZE*drand()); z = drand();}
if((gen9->n[m])%2 == 0){gen9->n[m+1] = (gen9->n[m]) + 1;} if((gen9->n[m])%2 == 1){gen9->n[m+1] = (gen9->n[m]) - 1;}
}
}

//Create next generation
j = gen9->n[m];if(gen9->np[m] < rec){if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0;while(position[l] <= gen9->brk[m]){selodd[m][l] = seleven[j][l]; ++l; if(l == N){break;}} memcpy(&selodd[m][l], &seleven[k][l], 2*(N-l));}
else{memcpy(&selodd[m][0],&seleven[j][0],2*N);}
j = gen9->n[m+1];if(gen9->np[m+1] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen9->brk[m+1]){selodd[m+1][l]=seleven[j][l];++l;if(l==N){break;}}memcpy(&selodd[m+1][l],&seleven[k][l],2*(N-l));}
else{memcpy(&selodd[m+1][0],&seleven[j][0],2*N);}

//Mutate at selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m][j]   == 0){selodd[m][j]   = 1;continue;}if(selodd[m][j]   == 1){selodd[m][j]   = 0;continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m+1][j] == 0){selodd[m+1][j] = 1;continue;}if(selodd[m+1][j] == 1){selodd[m+1][j] = 0;continue;}}

//Calculate fitness
ofit[m] = 1; ofit[m+1] = 1;
for(l=0; l < N; l++){
if(selodd[m][l]      != selodd[m+1][l]     ){if(FITNESS == 1){ofit[m] += H[l];  ofit[m+1] += H[l];}  else{ofit[m] *= H[l];  ofit[m+1]*= H[l];}}
if(selodd[m][l] == 1 && selodd[m+1][l] == 1){if(FITNESS == 1){ofit[m] += S[l];  ofit[m+1] += S[l];}  else{ofit[m] *= S[l];  ofit[m+1]*= S[l];}}
if(selodd[m][l] == 0 && selodd[m+1][l] == 0){if(FITNESS == 1){ofit[m] += NU[l]; ofit[m+1] += NU[l];} else{ofit[m] *= NU[l]; ofit[m+1]*= NU[l];}}
}

if(ofit[m] > maxfit){maxfit = ofit[m];}
m = m + 2;}

for(m=0;m<POPSIZE;m++){ofit[m] /= maxfit;}

memset(&n2count[0],'\0',4*POPSIZE);memset(&n2r[0],'\0',4*POPSIZE);
memset(&n3count[0],'\0',4*POPSIZE);memset(&n3r[0],'\0',4*POPSIZE);
memset(&n4count[0],'\0',4*POPSIZE);memset(&n4r[0],'\0',4*POPSIZE);
memset(&n5count[0],'\0',4*POPSIZE);memset(&n5r[0],'\0',4*POPSIZE);
memset(&n6count[0],'\0',4*POPSIZE);memset(&n6r[0],'\0',4*POPSIZE);
memset(&n7count[0],'\0',4*POPSIZE);memset(&n7r[0],'\0',4*POPSIZE);
memset(&n8count[0],'\0',4*POPSIZE);memset(&n8r[0],'\0',4*POPSIZE);
memset(&n9count[0],'\0',4*POPSIZE);memset(&n9r[0],'\0',4*POPSIZE);




for(i=1;i<=GEN;i++){
//Create odd generation
if(i%2 == 1){ 

//Initialize length for chromosome array
memset(&size[0],'\0',4*POPSIZE);

//Initialize the skipping array
memset(&indic[0],'\0',4*POPSIZE);

//Initialize recount array
memset(&recount[0],'\0',4*POPSIZE);


for(m=0;m<POPSIZE;m++){
temp1[m] = gen2->n[m];n2count[temp1[m]] = 1; 
if(gen2->np[m]<rec){
if(gen2->n[m]%2==0){n2r[gen2->n[m]+1] = 1;}
if(gen2->n[m]%2==1){n2r[gen2->n[m]-1] = 1;}
}current[m] = m;} 
  
for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen3->n[m]];n3count[temp2[m]] = 1; 
if(gen3->np[m]<rec){
if(gen3->n[m]%2==0){n3r[temp1[gen3->n[m]+1]] = 1;}
if(gen3->n[m]%2==1){n3r[temp1[gen3->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen4->n[m]];n4count[temp1[m]] = 1; 
if(gen4->np[m]<rec){
if(gen4->n[m]%2==0){n4r[temp2[gen4->n[m]+1]] = 1;}
if(gen4->n[m]%2==1){n4r[temp2[gen4->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen5->n[m]];n5count[temp2[m]] = 1; 
if(gen5->np[m]<rec){
if(gen5->n[m]%2==0){n5r[temp1[gen5->n[m]+1]] = 1;}
if(gen5->n[m]%2==1){n5r[temp1[gen5->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen6->n[m]];n6count[temp1[m]] = 1;
if(gen6->np[m]<rec){
if(gen6->n[m]%2==0){n6r[temp2[gen6->n[m]+1]] = 1;}
if(gen6->n[m]%2==1){n6r[temp2[gen6->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen7->n[m]];n7count[temp2[m]] = 1;
if(gen7->np[m]<rec){
if(gen7->n[m]%2==0){n7r[temp1[gen7->n[m]+1]] = 1;}
if(gen7->n[m]%2==1){n7r[temp1[gen7->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen8->n[m]];n8count[temp1[m]] = 1;
if(gen8->np[m]<rec){
if(gen8->n[m]%2==0){n8r[temp2[gen8->n[m]+1]] = 1;}
if(gen8->n[m]%2==1){n8r[temp2[gen8->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen9->n[m]];n9count[temp2[m]] = 1;
if(gen9->np[m]<rec){
if(gen9->n[m]%2==0){n9r[temp1[gen9->n[m]+1]] = 1;}
if(gen9->n[m]%2==1){n9r[temp1[gen9->n[m]-1]] = 1;}
}}


m=0;
while(m<POPSIZE){
flag1=0;
while(ZERO==0){
if(n2count[m]==0){flag1=1;break;}
if(n3r[m]==0){
if(n3count[m]==0){flag1=1;break;}
if(n4r[m]==0){
if(n4count[m]==0){flag1=1;break;}
if(n5r[m]==0){
if(n5count[m]==0){flag1=1;break;}
if(n6r[m]==0){
if(n6count[m]==0){flag1=1;break;}
if(n7r[m]==0){
if(n7count[m]==0){flag1=1;break;}
if(n8r[m]==0){
if(n8count[m]==0){flag1=1;break;}
if(n9r[m]==0){
if(n9count[m]==0){flag1=1;break;}
}}}}}}}
break;}

flag2=0;
while(ZERO==0){
if(n2count[m+1]==0){flag2=1;break;}
if(n3r[m+1]==0){
if(n3count[m+1]==0){flag2=1;break;}
if(n4r[m+1]==0){
if(n4count[m+1]==0){flag2=1;break;}
if(n5r[m+1]==0){
if(n5count[m+1]==0){flag2=1;break;}
if(n6r[m+1]==0){
if(n6count[m+1]==0){flag2=1;break;}
if(n7r[m+1]==0){
if(n7count[m+1]==0){flag2=1;break;}
if(n8r[m+1]==0){
if(n8count[m+1]==0){flag2=1;break;}
if(n9r[m+1]==0){
if(n9count[m+1]==0){flag2=1;break;}
}}}}}}}
break;}

if(flag1 == 1 && flag2    == 1){indic[m]=1;indic[m+1]=1;}
else{
if(flag1 == 1 && n2r[m]   == 0){indic[m]=1;}
if(flag2 == 1 && n2r[m+1] == 0){indic[m+1]=1;}
}

if(indic[m] == 0   || (i%DELETE >= DELETE-10 || i%DELETE == 0 || i >= GEN-10)){
p = gen1->np[m];j = gen1->n[m];
if(p < rec){if((j&1)==0){k = j+1;}else{k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&odd[m].chr->q[0],&odd[current[j]-POPSIZE].chr->q[0],4*siz[j]);size[m]=siz[j];}
if(current[j]==j){temp = odd[m].chr;odd[m].chr = even[j].chr;even[j].chr = temp;current[j]=m + POPSIZE;size[m] = siz[j];}
}}
  
if(indic[m+1] == 0 || (i%DELETE >= DELETE-10 || i%DELETE == 0 || i >= GEN-10)){
p = gen1->np[m+1];j = gen1->n[m+1];
if(p < rec){if((j&1)==0){k = j+1;}else{k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&odd[m+1].chr->q[0],&odd[current[j]-POPSIZE].chr->q[0],4*siz[j]);size[m+1]=siz[j];}
if(current[j]==j){temp = odd[m+1].chr;odd[m+1].chr = even[j].chr;even[j].chr = temp;current[j]=m+1+POPSIZE;size[m+1] = siz[j];}
}}
m=m+2;}



for(m=0;m<POPSIZE;m++){
p = gen1->np[m];j = gen1->n[m];
if(i%DELETE > 0  && i%DELETE < DELETE-10 && i < GEN-10){if(indic[m]==1){continue;}}

if(p<rec){
r = gen1->brk[m];if((j&1)==0){k = j+1;}else{k = j-1;}--recount[(j+k)/2];
if(current[j]==j && siz[j] > 0){
if(r>=even[j].chr->q[siz[j]-1]){
if(recount[(j+k)/2]==0){temp=odd[m].chr;odd[m].chr=even[j].chr;even[j].chr=temp;size[m]=siz[j];}
else{memcpy(&odd[m].chr->q[0],&even[j].chr->q[0],4*siz[j]);size[m]=siz[j];}
}
  
else{
start=0;end=siz[j]-1;trk=end;ind=0;
while(trk > 0 && trk <=siz[j]-1){
if(r >=even[j].chr->q[trk-1] && r <even[j].chr->q[trk]){ind=1;break;}
if(r >=even[j].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[j].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){
if(recount[(j+k)/2]==0){temp=odd[m].chr;odd[m].chr=even[j].chr;even[j].chr=temp;size[m]=trk;}
else{memcpy(&odd[m].chr->q[0],&even[j].chr->q[0],4*trk);size[m]=trk;}
}}
}
 
if(current[j]!=j && siz[j] > 0){
trp=current[j]-POPSIZE;  
if(r>=odd[trp].chr->q[siz[j]-1]){memcpy(&odd[m].chr->q[0],&odd[trp].chr->q[0],4*siz[j]);size[m]=siz[j];}
else{
start=0;end=siz[j]-1;trk=end;ind=0;
while(trk > 0 && trk<=siz[j]-1){
if(r >=odd[trp].chr->q[trk-1] && r <odd[trp].chr->q[trk]){ind=1;break;}
if(r >=odd[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&odd[m].chr->q[0],&odd[trp].chr->q[0],4*trk);size[m]=trk;}}}

if(current[k]==k && siz[k]>0){
if(r<even[k].chr->q[0]){memcpy(&odd[m].chr->q[size[m]],&even[k].chr->q[0],4*siz[k]);size[m]=size[m]+siz[k];}
else{
start=0;end=siz[k]-1;trk=end;ind=0;
while(trk > 0 && trk <=siz[k]-1){
if(r >=even[k].chr->q[trk-1] && r <even[k].chr->q[trk]){ind=1;break;}
if(r >=even[k].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[k].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&odd[m].chr->q[size[m]],&even[k].chr->q[trk],4*(siz[k]-trk));size[m]=size[m]+siz[k]-trk;}}}

if(current[k]!=k && siz[k]>0){trp=current[k]-POPSIZE;
if(r<odd[trp].chr->q[0]){memcpy(&odd[m].chr->q[size[m]],&odd[trp].chr->q[0],4*siz[k]);size[m]=size[m]+siz[k];}
else{
start=0;end=siz[k]-1;trk=end;ind=0;
while(trk > 0 && trk<=siz[k]-1){
if(r >=odd[trp].chr->q[trk-1] && r <odd[trp].chr->q[trk]){ind=1;break;}
if(r >=odd[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&odd[m].chr->q[size[m]],&odd[trp].chr->q[trk],4*(siz[k]-trk));size[m]=size[m]+siz[k]-trk;}}}
}}  




//Insert new mutations into the chromosome arrays
for(m=0; m < POPSIZE; m++){ 

if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m] == 1){continue;}}

//INSERT NEW MUTATIONS
k = poisson(mu);
if(k + size[m] > array_size){cout<<"Out of Memory. Increasing array size and rerunning.\n";break;}

for(l=0;l<k;l++){
trp = int(SEQLENGTH*drand());while(multhit[trp] == 1 || selected[trp] == 1){trp = int(SEQLENGTH*drand());} multhit[trp] = 1;
if((size[m] > 0 && trp>=odd[m].chr->q[size[m]-1]) || size[m]==0){odd[m].chr->q[size[m]]=trp;}
else{
if(trp <= odd[m].chr->q[0]){memmove(&odd[m].chr->q[1],&odd[m].chr->q[0],4*size[m]);odd[m].chr->q[0]=trp;}
else{
start = 0; end = size[m]-1;trk = end; ind = 0;
while(trk > 0 && trk<=size[m]-1){
if(trp > odd[m].chr->q[trk-1] && trp <=odd[m].chr->q[trk]){ind=1;break;} 
if(trp > odd[m].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;} 
if(trp <= odd[m].chr->q[trk]){end=trk;trk=(start+end)/2;continue;}} 
if(ind == 1){memmove(&odd[m].chr->q[trk+1],&odd[m].chr->q[trk],4*(size[m]-trk));odd[m].chr->q[trk]=trp;}
if(ind == 0){cout<<"Error\n";return(0);}}
}++size[m];}
 
}

if(m < POPSIZE){flag = 1;break;}


//Update the future genealogy by one more generation
gtemp = gen1; gen1 = gen2; gen2 = gen3; gen3 = gen4; gen4 = gen5; gen5 = gen6; gen6 = gen7; gen7 = gen8; gen8 = gen9; gen9 = gtemp;

m = 0; maxfit = 0;
while(m < POPSIZE){
//Choose chromosomes according to fitness
gen9->np[m] = drand();gen9->np[m+1] = drand(); x = drand();
if(gen9->np[m] < rec){gen9->brk[m] = int((SEQLENGTH-1)*drand());}if(gen9->np[m+1] < rec){gen9->brk[m+1] = int((SEQLENGTH-1)*drand());}

//If not created by selfing
if(x >= s){
gen9->n[m]   = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen9->n[m]]  ){gen9->n[m]   = int(POPSIZE*drand()); z = drand();}
gen9->n[m+1] = int(POPSIZE*drand()); z = drand(); while(z >= ofit[gen9->n[m+1]]){gen9->n[m+1] = int(POPSIZE*drand()); z = drand();}
}

//If created by selfing
else{
y = drand();
if(y <= 0.50){
gen9->n[m]   = int(POPSIZE*drand()); z = drand();while(z >= ofit[gen9->n[m]]){gen9->n[m]=int(POPSIZE*drand());z=drand();}
gen9->n[m+1] = gen9->n[m];}else{gen9->n[m]=int(POPSIZE*drand());z=drand();while(z >= ofit[gen9->n[m]]){gen9->n[m]=int(POPSIZE*drand());z=drand();}
if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}if((gen9->n[m])%2==1){gen9->n[m+1] = (gen9->n[m]) - 1;}}
}

//Create next generation
j = gen9->n[m];if(gen9->np[m] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen9->brk[m]){seleven[m][l]=selodd[j][l];++l;if(l==N){break;}}memcpy(&seleven[m][l],&selodd[k][l],2*(N-l));}
else{memcpy(&seleven[m][0],&selodd[j][0],2*N);}
j = gen9->n[m+1];if(gen9->np[m+1] < rec){if((j&1)==0){k = j+1;}else{k = j-1;}
l=0;while(position[l] <= gen9->brk[m+1]){seleven[m+1][l]=selodd[j][l];++l;if(l==N){break;}}memcpy(&seleven[m+1][l],&selodd[k][l],2*(N-l));}
else{memcpy(&seleven[m+1][0],&selodd[j][0],2*N);}

//Mutate
k = poisson(V);
for(l=0;l<k;l++){j = int(N*drand()); if(seleven[m][j]   == 0){seleven[m  ][j] = 1;continue;} if(seleven[m  ][j] == 1){seleven[m][j]   = 0;continue;}}
k = poisson(V);
for(l=0;l<k;l++){j = int(N*drand()); if(seleven[m+1][j] == 0){seleven[m+1][j] = 1;continue;} if(seleven[m+1][j] == 1){seleven[m+1][j] = 0;continue;}}

//Calculate fitness
efit[m] = 1; efit[m+1] = 1;
for(l=0;l<N;l++){
if(seleven[m][l]      !=  seleven[m+1][l]     ){if(FITNESS == 1){efit[m] += H[l]; efit[m+1] += H[l];}  else{efit[m] *= H[l];  efit[m+1] *= H[l]; }}
if(seleven[m][l] == 1 &&  seleven[m+1][l] == 1){if(FITNESS == 1){efit[m] += S[l]; efit[m+1] += S[l];}  else{efit[m] *= S[l];  efit[m+1] *= S[l]; }}
if(seleven[m][l] == 0 &&  seleven[m+1][l] == 0){if(FITNESS == 1){efit[m] += NU[l];efit[m+1] += NU[l];} else{efit[m] *= NU[l]; efit[m+1] *= NU[l];}}
}

if(efit[m] > maxfit){maxfit = efit[m];}
m = m + 2;}

//Scale fitness in terms of maximum fitness
for(m=0;m<POPSIZE;m++){efit[m] /= maxfit;}

memset(&n2count[0],'\0',4*POPSIZE);memset(&n2r[0],'\0',4*POPSIZE);
memset(&n3count[0],'\0',4*POPSIZE);memset(&n3r[0],'\0',4*POPSIZE);
memset(&n4count[0],'\0',4*POPSIZE);memset(&n4r[0],'\0',4*POPSIZE);
memset(&n5count[0],'\0',4*POPSIZE);memset(&n5r[0],'\0',4*POPSIZE);
memset(&n6count[0],'\0',4*POPSIZE);memset(&n6r[0],'\0',4*POPSIZE);
memset(&n7count[0],'\0',4*POPSIZE);memset(&n7r[0],'\0',4*POPSIZE);
memset(&n8count[0],'\0',4*POPSIZE);memset(&n8r[0],'\0',4*POPSIZE);
memset(&n9count[0],'\0',4*POPSIZE);memset(&n9r[0],'\0',4*POPSIZE);
}



 
if(i%2 == 0){ 
//Initialize length of chromosome array
memset(&siz[0],'\0',4*POPSIZE);

//Initialize skipping array
memset(&indic[0],'\0',4*POPSIZE);

//Initialize recount array
memset(&recount[0],'\0',4*POPSIZE);


for(m=0;m<POPSIZE;m++){
temp1[m] = gen2->n[m];n2count[temp1[m]] = 1; 
if(gen2->np[m] < rec){
if(gen2->n[m]%2 == 0){n2r[gen2->n[m]+1] = 1;}
if(gen2->n[m]%2 == 1){n2r[gen2->n[m]-1] = 1;}
} current[m] = m;} 

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen3->n[m]];n3count[temp2[m]] = 1; 
if(gen3->np[m]<rec){
if(gen3->n[m]%2==0){n3r[temp1[gen3->n[m]+1]] = 1;}
if(gen3->n[m]%2==1){n3r[temp1[gen3->n[m]-1]] = 1;}
}} 
  
for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen4->n[m]];n4count[temp1[m]] = 1; 
if(gen4->np[m]<rec){
if(gen4->n[m]%2==0){n4r[temp2[gen4->n[m]+1]] = 1;}
if(gen4->n[m]%2==1){n4r[temp2[gen4->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen5->n[m]];n5count[temp2[m]] = 1; 
if(gen5->np[m]<rec){
if(gen5->n[m]%2==0){n5r[temp1[gen5->n[m]+1]] = 1;}
if(gen5->n[m]%2==1){n5r[temp1[gen5->n[m]-1]] = 1;}
}} 

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen6->n[m]];n6count[temp1[m]] = 1;
if(gen6->np[m]<rec){
if(gen6->n[m]%2==0){n6r[temp2[gen6->n[m]+1]] = 1;}
if(gen6->n[m]%2==1){n6r[temp2[gen6->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen7->n[m]];n7count[temp2[m]] = 1;
if(gen7->np[m]<rec){
if(gen7->n[m]%2==0){n7r[temp1[gen7->n[m]+1]] = 1;}
if(gen7->n[m]%2==1){n7r[temp1[gen7->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp1[m] = temp2[gen8->n[m]];n8count[temp1[m]] = 1;
if(gen8->np[m]<rec){
if(gen8->n[m]%2==0){n8r[temp2[gen8->n[m]+1]] = 1;}
if(gen8->n[m]%2==1){n8r[temp2[gen8->n[m]-1]] = 1;}
}}

for(m=0;m<POPSIZE;m++){
temp2[m] = temp1[gen9->n[m]];n9count[temp2[m]] = 1;
if(gen9->np[m]<rec){
if(gen9->n[m]%2==0){n9r[temp1[gen9->n[m]+1]] = 1;}
if(gen9->n[m]%2==1){n9r[temp1[gen9->n[m]-1]] = 1;}
}}




m=0;
while(m < POPSIZE){
flag1 = 0;
while(ZERO==0){
if(n2count[m]==0){flag1=1;break;}
if(n3r[m]==0){
if(n3count[m]==0){flag1=1;break;}
if(n4r[m]==0){
if(n4count[m]==0){flag1=1;break;}
if(n5r[m]==0){
if(n5count[m]==0){flag1=1;break;}
if(n6r[m]==0){
if(n6count[m]==0){flag1=1;break;}
if(n7r[m]==0){
if(n7count[m]==0){flag1=1;break;}
if(n8r[m]==0){
if(n8count[m]==0){flag1=1;break;}
if(n9r[m]==0){
if(n9count[m]==0){flag1=1;break;}
}}}}}}}
break;}


flag2=0;
while(ZERO==0){
if(n2count[m+1]==0){flag2=1;break;}
if(n3r[m+1]==0){
if(n3count[m+1]==0){flag2=1;break;}
if(n4r[m+1]==0){
if(n4count[m+1]==0){flag2=1;break;}
if(n5r[m+1]==0){
if(n5count[m+1]==0){flag2=1;break;}
if(n6r[m+1]==0){
if(n6count[m+1]==0){flag2=1;break;}
if(n7r[m+1]==0){
if(n7count[m+1]==0){flag2=1;break;}
if(n8r[m+1]==0){
if(n8count[m+1]==0){flag2=1;break;}
if(n9r[m+1]==0){
if(n9count[m+1]==0){flag2=1;break;}
}}}}}}}
break;}

if(flag1==1 && flag2==1){indic[m]=1;indic[m+1]=1;}
else{
if(flag1==1 && n2r[m]==0){indic[m]=1;}
if(flag2==1 && n2r[m+1]==0){indic[m+1]=1;}}
  
if(indic[m]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
p = gen1->np[m];j = gen1->n[m];if(p < rec){if((j&1)==0){k = j+1;}else{k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&even[m].chr->q[0],&even[current[j]-POPSIZE].chr->q[0],4*size[j]);siz[m]=size[j];}
if(current[j]==j){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+POPSIZE;siz[m] = size[j];}
}}

if(indic[m+1]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
p = gen1->np[m+1];j = gen1->n[m+1];if(p < rec){if((j&1)==0){k = j+1;}else{k = j-1;}++recount[(j+k)/2];}
else{
if(current[j]!=j){memcpy(&even[m+1].chr->q[0],&even[current[j]-POPSIZE].chr->q[0],4*size[j]);siz[m+1]=size[j];}
if(current[j]==j){temp = even[m+1].chr;even[m+1].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+1+POPSIZE;siz[m+1] = size[j];}
}}
m=m+2;}

for(m=0;m<POPSIZE;m++){
p = gen1->np[m];j = gen1->n[m];
if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}}

if(p<rec){
r = gen1->brk[m];if((j&1)==0){k = j+1;}else{k = j-1;}--recount[(j+k)/2];  
if(current[j]==j && size[j]>0){
if(r>=odd[j].chr->q[size[j]-1]){
if(recount[(j+k)/2]==0){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j] = m+POPSIZE;siz[m] = size[j];}
else{memcpy(&even[m].chr->q[0],&odd[j].chr->q[0],4*size[j]);siz[m]=size[j];}
}
  
else{
start=0;end=size[j]-1;trk=end;ind=0;
while(trk > 0 && trk <=size[j]-1){
if(r >=odd[j].chr->q[trk-1] && r <odd[j].chr->q[trk]){ind=1;break;}
if(r >=odd[j].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[j].chr->q[trk]){end=trk;trk=(start+end)/2;}
}

if(ind==1){
if(recount[(j+k)/2]==0){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+POPSIZE;siz[m] = trk;}  
else{memcpy(&even[m].chr->q[0],&odd[j].chr->q[0],4*trk);siz[m]=trk;}
}}
}

if(current[j]!=j && size[j]>0){
trp=current[j]-POPSIZE;
if(r>=even[trp].chr->q[size[j]-1]){memcpy(&even[m].chr->q[0],&even[trp].chr->q[0],4*size[j]);siz[m]=size[j];}
else{
start=0;end=size[j]-1;trk=end;ind=0;
while(trk > 0 && trk<=size[j]-1){
if(r >=even[trp].chr->q[trk-1] && r <even[trp].chr->q[trk]){ind=1;break;}
if(r >=even[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&even[m].chr->q[0],&even[trp].chr->q[0],4*trk);siz[m]=trk;}}}

if(current[k]==k && size[k]>0){
if(r<odd[k].chr->q[0]){memcpy(&even[m].chr->q[siz[m]],&odd[k].chr->q[0],4*size[k]);siz[m]=siz[m]+size[k];}
else{
start=0;end=size[k]-1;trk=end;ind=0;
while(trk > 0 && trk <=size[k]-1){
if(r >=odd[k].chr->q[trk-1] && r <odd[k].chr->q[trk]){ind=1;break;}
if(r >=odd[k].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <odd[k].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&even[m].chr->q[siz[m]],&odd[k].chr->q[trk],4*(size[k]-trk));siz[m]=siz[m]+size[k]-trk;}}}

if(current[k]!=k && size[k]>0){
trp=current[k]-POPSIZE;
if(r<even[trp].chr->q[0]){memcpy(&even[m].chr->q[siz[m]],&even[trp].chr->q[0],4*size[k]);siz[m]=siz[m]+size[k];}
else{
start=0;end=size[k]-1;trk=end;ind=0;
while(trk > 0 && trk<=size[k]-1){
if(r >=even[trp].chr->q[trk-1] && r <even[trp].chr->q[trk]){ind=1;break;}
if(r >=even[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(r <even[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}}
if(ind==1){memcpy(&even[m].chr->q[siz[m]],&even[trp].chr->q[trk],4*(size[k]-trk));siz[m]=siz[m]+size[k]-trk;}}}
}}
  
  
  
//INSERT NEW MUTATIONS    
for(m=0;m<POPSIZE;m++){

if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m] == 1){continue;}}

k = poisson(mu);
if(k + siz[m] > array_size){cout<<"Out of Memory. Increasing array size and rerunning.\n"; break;}

for(l=0;l<k;l++){
trp = int(SEQLENGTH*drand());while(multhit[trp]  ==  1 || selected[trp]  == 1){trp = int(SEQLENGTH*drand());}multhit[trp] = 1;
if((siz[m]> 0 && trp>=even[m].chr->q[siz[m]-1]) || siz[m]==0){even[m].chr->q[siz[m]]=trp;}
else{
if(trp<=even[m].chr->q[0]){memmove(&even[m].chr->q[1],&even[m].chr->q[0],4*siz[m]);even[m].chr->q[0]=trp;}
else{
start=0;end=siz[m]-1;trk=end;ind=0;
while(trk > 0 && trk<=siz[m]-1){
if(trp > even[m].chr->q[trk-1] && trp <=even[m].chr->q[trk]){ind=1;break;}
if(trp > even[m].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}continue;}
if(trp <=even[m].chr->q[trk]){end=trk;trk=(start+end)/2;continue;}
}

if(ind  == 1){memmove(&even[m].chr->q[trk+1],&even[m].chr->q[trk],4*(siz[m]-trk));even[m].chr->q[trk]=trp;}
if(ind  == 0){cout<<"Error\n";return(0);}}
}++siz[m];}
}

if(m < POPSIZE){flag = 1; break;}

//Update the future genealogy by one more generation
gtemp = gen1; gen1 = gen2; gen2 = gen3; gen3 = gen4; gen4 = gen5; gen5 = gen6; gen6 = gen7; gen7 = gen8; gen8 = gen9; gen9 = gtemp;

m = 0; maxfit = 0;
while(m < POPSIZE){
//Choose chromosomes according to fitness
gen9->np[m] = drand(); gen9->np[m+1] = drand(); 
if(gen9->np[m] < rec){gen9->brk[m] = int((SEQLENGTH-1)*drand());}if(gen9->np[m+1] < rec){gen9->brk[m+1] = int((SEQLENGTH-1)*drand());}

x = drand();
//If not created by selfing
if(x >= s){
gen9->n[m]   = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen9->n[m]]  ){gen9->n[m]   = int(POPSIZE*drand()); z = drand();}
gen9->n[m+1] = int(POPSIZE*drand()); z = drand(); while(z >= efit[gen9->n[m+1]]){gen9->n[m+1] = int(POPSIZE*drand()); z = drand();}
}

//If created by selfing
else{
y = drand();
if(y <=0.50){
gen9->n[m] = int(POPSIZE*drand());z = drand(); while(z >= efit[gen9->n[m]]){gen9->n[m]=int(POPSIZE*drand());z = drand();}
gen9->n[m+1] = gen9->n[m];
}
else{gen9->n[m]=int(POPSIZE*drand());z=drand();while(z >= efit[gen9->n[m]]){gen9->n[m]=int(POPSIZE*drand());z = drand();}
if((gen9->n[m])%2 == 0){gen9->n[m+1] = (gen9->n[m])+1;}if((gen9->n[m])%2 == 1){gen9->n[m+1] = (gen9->n[m]) - 1;}}
}

//Create next generation 
j = gen9->n[m]; 
if(gen9->np[m] < rec){
if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0; while(position[l] <= gen9->brk[m]  ){selodd[m][l]   =  seleven[j][l]; ++l; if(l == N){break;}} memcpy(&selodd[m][l],&seleven[k][l],2*(N-l));
}
else{memcpy(&selodd[m][0],   &seleven[j][0], 2*N);}

j = gen9->n[m+1];
if(gen9->np[m+1] < rec){
if((j&1) == 0){k = j + 1;}else{k = j - 1;}
l = 0; while(position[l] <= gen9->brk[m+1]){selodd[m+1][l]  =  seleven[j][l]; ++l; if(l == N){break;}} memcpy(&selodd[m+1][l],&seleven[k][l],2*(N-l));
}
else{memcpy(&selodd[m+1][0], &seleven[j][0], 2*N);}

//Mutate the selected sites
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m][j]   == 0){selodd[m][j] = 1;   continue;}if(selodd[m][j]   == 1){selodd[m][j] = 0;  continue;}}
k = poisson(V);for(l=0;l<k;l++){j = int(N*drand());if(selodd[m+1][j] == 0){selodd[m+1][j] = 1; continue;}if(selodd[m+1][j] == 1){selodd[m+1][j] = 0;continue;}}

//Calculate fitness
ofit[m] = 1; ofit[m+1] = 1;
for(l=0;l<N;l++){
if(selodd[m][l]      != selodd[m+1][l]     ){if(FITNESS==1){ofit[m]+= H[l];  ofit[m+1]+=  H[l];} else{ofit[m] *= H[l];  ofit[m+1] *=  H[l];}}
if(selodd[m][l] == 1 && selodd[m+1][l] == 1){if(FITNESS==1){ofit[m]+= S[l];  ofit[m+1]+=  S[l];} else{ofit[m] *= S[l];  ofit[m+1] *=  S[l];}}
if(selodd[m][l] == 0 && selodd[m+1][l] == 0){if(FITNESS==1){ofit[m]+= NU[l]; ofit[m+1]+= NU[l];} else{ofit[m] *= NU[l]; ofit[m+1] *= NU[l];}}
}

if(ofit[m] > maxfit){maxfit = ofit[m];}
m = m + 2;}

//Scale fitness in terms of maximum fitness
for(m=0; m<POPSIZE; m++){ofit[m] /= maxfit;}


memset(&n2count[0],'\0',4*POPSIZE);memset(&n2r[0],'\0',4*POPSIZE);
memset(&n3count[0],'\0',4*POPSIZE);memset(&n3r[0],'\0',4*POPSIZE);
memset(&n4count[0],'\0',4*POPSIZE);memset(&n4r[0],'\0',4*POPSIZE);
memset(&n5count[0],'\0',4*POPSIZE);memset(&n5r[0],'\0',4*POPSIZE);
memset(&n6count[0],'\0',4*POPSIZE);memset(&n6r[0],'\0',4*POPSIZE);
memset(&n7count[0],'\0',4*POPSIZE);memset(&n7r[0],'\0',4*POPSIZE);
memset(&n8count[0],'\0',4*POPSIZE);memset(&n8r[0],'\0',4*POPSIZE);
memset(&n9count[0],'\0',4*POPSIZE);memset(&n9r[0],'\0',4*POPSIZE);
}




//REMOVE FIXED MUTATIONS FROM THE CHROMOSOME ARRAYS
if(i%DELETE == 0 || i == GEN){

if(i%2 ==  0){
memset(&indicator[0],'\0',4*SEQLENGTH);
for(m=0;m<POPSIZE;m++){
if(siz[m] > array_size){cout<<"Out of Memory. Increase array size of q.\n"; return(0);}
for(l=0; l<siz[m]; l++){++indicator[even[m].chr->q[l]];}
}

for(l=0; l < SEQLENGTH; l++){
if(indicator[l] == 0 || indicator[l] == POPSIZE){multhit[l] = 0;}
if(indicator[l] > POPSIZE || indicator[l] < 0){cout<<i<<" "<<"error\n"; return(0);}
}

//Remove fixed locations
int count = 0;
for(m=0;m<POPSIZE;m++){
k = 0; trp = siz[m];l = 0;
while(l < trp){
trk = even[m].chr->q[l];
if(indicator[trk] == POPSIZE){++l; continue;}
else{even[m].chr->q[k] = trk; ++k; ++l; continue;}
} 
siz[m] = k; count = count + siz[m];}
}

else{
memset(&indicator[0],'\0',4*SEQLENGTH);
for(m=0; m<POPSIZE; m++){
if(size[m] > array_size){cout<<"Out of Memory. Increase array size of q.\n"; return(0);}
for(l=0;l<size[m];l++){++indicator[odd[m].chr->q[l]];}
}

for(l=0;l<SEQLENGTH; l++){
if(indicator[l] == 0 || indicator[l] == POPSIZE){multhit[l] = 0;}
if(indicator[l] > POPSIZE || indicator[l] < 0){cout<<i<<" "<<"error\n"; return(0);}
}

int count = 0;
for(m=0; m<POPSIZE; m++){
k = 0;trp = size[m]; l = 0;
while(l< trp){
trk = odd[m].chr->q[l];
if(indicator[trk] == POPSIZE){++l;       continue;}
else{odd[m].chr->q[k] = trk; ++k; ++l;   continue;}
} size[m] = k; count = count + size[m];}
} 
}


//Copy the alleles of selected sites in the final (GEN th) generation
if(i == GEN - 8){
if(i%2 == 0){memcpy(&fsample[0][0], &seleven[0][0], 2*POPSIZE*N);}
else        {memcpy(&fsample[0][0], &selodd[0][0],  2*POPSIZE*N);}
}

}




//If number of mutations exceed arrays_size
if(flag == 1){
for(m = 0; m<POPSIZE; m++){
delete odd[m].chr->q;
delete odd[m].chr;
delete even[m].chr->q;
delete even[m].chr;
}
continue;}




//Get the locations of polymorphic sites
memset(&indicator[0],'\0',4*SEQLENGTH);
if(GEN%2 == 0){for(m=0;m<SSIZE;m++){for(i=0;i<siz[2*m] ; i++){++indicator[even[2*m].chr->q[i]];}}}
else          {for(m=0;m<SSIZE;m++){for(i=0;i<size[2*m]; i++){++indicator[odd[2*m].chr->q[i]]; }}}

tmpp = new int[SEQLENGTH]; memset(&tmpp[0],'\0',4*SEQLENGTH); 

mut = 0;
for(i=0;i<SEQLENGTH;i++){
if(selected[i] == 1){
for(k=0;k<N;k++){if(position[k] == i){break;}}
assert(k < N);

l = 0; for(j=0; j<SSIZE; j++){l += fsample[2*j][k];}
if(l >= 1 && l <= SSIZE-1){tmpp[i] = 1; ++mut;}
}

else{
if(indicator[i] >= 1 && indicator[i] <= SSIZE-1){tmpp[i] = 1; ++mut;}
}
}


//Output polymorphic locations for the sample in the final generation
cout<<"Sample size: "<<SSIZE<<"\n";
cout<<"SNP Count: "<<mut<<" "<<"\n";
cout<<"Positions: ";for(i=0;i<SEQLENGTH;i++){if(tmpp[i] == 1){cout<<i+1<<" ";}}cout<<"\n";

for(j=0;j<SSIZE;j++){
for(i=0;i<SEQLENGTH;i++){
if(tmpp[i] == 0){continue;}

if(selected[i] == 1){
for(k=0;k<N;k++){if(position[k] == i){break;}}
if(k<N){cout<<fsample[2*j][k];}
else{cout<<"ERR\n";return(0);}
}

else{
if(GEN%2 == 0){int flag = 0; for(l=0; l<siz[2*j]; l++){if(even[2*j].chr->q[l] == i){cout<<"1";flag = 1;break;}} if(flag == 0){cout<<"0";}}
else{int flag = 0; for(l=0; l<size[2*j]; l++){if(odd[2*j].chr->q[l] == i){cout<<"1"; flag = 1; break;}} if(flag == 0){cout<<"0";}}
}

}cout<<"\n";}
cout<<"\n";


delete tmpp;


//Delete chromosome arrays
for(m=0;m<POPSIZE;m++){
delete odd[m].chr->q;
delete odd[m].chr;
delete even[m].chr->q;
delete even[m].chr;
}

finished = 1;}
  


}

}



