/*
 * RandomizedMotifSearch
 * Delete the latest "\n" from input file
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DEBUG

typedef unsigned char    UINT8;
typedef unsigned int     UINT32;


#define MAX_ARRAY_SIZE   1000000


#define SCAN_KMER        0
#define SCAN_NUMOFDNA    1
#define SCAN_DNA         2


typedef struct {
    UINT8 *src;
    UINT32 size_of_each_dna;
    UINT32 kmer;
    UINT32 num_of_dna;
    UINT32 total_size;
}DNA_INFO_t;

UINT32 pMostProbableKmer[2][100];//Store the most probability position.
UINT32 alpha[256] = {0};


/*
 * Read and parse data.
 */
bool read_test_data(const char* file_name,DNA_INFO_t *dna_info,int maxsize) 
{    
    UINT32 phase = SCAN_KMER,count=0;
    UINT8 *src = dna_info->src;
    UINT8 c;
    FILE* file = fopen (file_name, "r");
    
    if (!file) return FALSE;
    if (!src)  return FALSE;
    
    dna_info->total_size = dna_info->num_of_dna = dna_info->kmer = 0;
    
    do {
        if(SCAN_DNA == phase){
            c = fgetc(file);
            if( (c == 'G') || (c == 'C') || (c=='A') || (c=='T')){
                src[dna_info->total_size ++] = c;                
            }else if((c == 0xd) || (c == 0xa)){
                continue;
            }
        }

        if(SCAN_KMER == phase){
            c = fgetc(file);
            if(c == ' '){
                phase = SCAN_NUMOFDNA;
                continue;
            }
            dna_info->kmer *= 10;
            dna_info->kmer += c-'0';
        }
        
        if(SCAN_NUMOFDNA == phase){
            c = fgetc(file);                     
            if((c == 0xd) || (c == 0xa)){
                if(0xa == c)
                    phase = SCAN_DNA;
                continue;
            }

            dna_info->num_of_dna *= 10;
            dna_info->num_of_dna += c-'0';
        }        

    } while (!feof (file) && count < maxsize);
    
    dna_info->size_of_each_dna = dna_info->total_size/dna_info->num_of_dna;
    
    fclose (file);
    
    return TRUE;
}



/*
 * k: KMER
 */
double caculate_pr(UINT8 *dna,UINT32 k)
{
    UINT32 i;
    double multiple = 1;
    double j;
    
    for(i=0,multiple=1;i<k;i++){
        
        j = PROFILE_PATTERN[i]+PROFILE_PATTERN[i+KMER_LENGTH]+PROFILE_PATTERN[i+KMER_LENGTH*2]+PROFILE_PATTERN[i+KMER_LENGTH*3];
        j = PROFILE_PATTERN[alpha[dna[i]]*KMER_LENGTH+i]/j;
                    
        multiple*=j;
    }

    return multiple;
}


/*
 * Score the matrix motifs.
 * dna_src: Source of DNA.
 * kmer:  KMER length.
 * num_dnas: Number of DNA.
 * size_dna: Size of single DNA.
 */
UINT32 score_mmotifs(DNA_INFO_t *dna_info)
{
    UINT32 i,j,k,l;
    UINT32 sum_A,sum_C,sum_G,sum_T,max_NU=0,sum_column=0,sum=0;
    UINT8 *dna_src = dna_info->src;
    UINT32 kmer = dna_info->kmer;
    UINT32 num_dnas = dna_info->num_of_dna;
    UINT32 size_dna = dna_info->size_of_each_dna;

    for(i=0;i<kmer;i++){
        sum_A = sum_C = sum_G = sum_T = 0;
        for(j=0;j<num_dnas;j++){
            k=pMostProbableKmer[0][j];
            if(dna_src[j*size_dna+k+i] == 'A') sum_A++;
            if(dna_src[j*size_dna+k+i] == 'C') sum_C++;
            if(dna_src[j*size_dna+k+i] == 'G') sum_G++;
            if(dna_src[j*size_dna+k+i] == 'T') sum_T++;
        }

        sum_column = sum_A+sum_C+sum_G+sum_T;

        l = (sum_A>sum_C)?sum_A:sum_C;
        l = (sum_G>l)?sum_G:l;
        l = (sum_T>l)?sum_T:l;
        sum_column -= l;
        sum+=sum_column;
    }
    return sum;
}

/*
 * Limit to 1000 times
 * Dna: DNA strings.
 * kmer: KMER.
 * num_dnas: Number of DNA strings.
 * size_dna: Size of DNA.
 */
void RandomizedMotifSearch(UINT8 *Dna,const UINT32 kmer,const UINT32 num_dnas,const UINT32 size_dna)
{
    UINT32 i,pseudo_cnt=0;
    UINT32 score=0;
    UINT32 i=0,j;
    double m_prob,c_prob;

    do{
        for(j=1;j<t;j++){        
            for(i=0,m_prob=0,pMostProbableKmer[0][j] = 0;i<=(size_dna-kmer);i++){              
                  c_prob = caculate_pr(Dna+i+j*size_dna,kmer);
                  if(c_prob>m_prob){
                      m_prob = c_prob;
                      //Update the max probabity position
                      pMostProbableKmer[0][j] = i;
                  }
            }
        }

        //Score
        if()

        pseudo_cnt++;
    }while(pseudo_cnt<1000)
}

void main()
{
    UINT32 i=0,j=0;
    UINT32 a_cnt,c_cnt,g_cnt,t_cnt;
    UINT32 pr_multi,max_pr=0;
    UINT32 best_score;
    UINT32 *profile_pattern;
    DNA_INFO_t dna_info;
    
    dna_info.src = malloc(MAX_ARRAY_SIZE);
    if(!read_test_data("dataset_158_9.txt",&dna_info,MAX_ARRAY_SIZE))
        return;
#ifdef DEBUG
    printf("num_d:%d dna_total:%d dna_sz:%d kmer:%d \n",NUMBER_OF_DNA,DNA_TOTAL_LENGTH,SIZE_OF_DNA,KMER_LENGTH);
#endif
    // Init hash index.
    alpha['A'] = 0;
    alpha['C'] = 1;
    alpha['G'] = 2;
    alpha['T'] = 3;
  
    //Init profile table
    PROFILE_PATTERN = malloc(sizeof(UINT32)*KMER_LENGTH*NUMBER_OF_DNA);
    memset(PROFILE_PATTERN,0,sizeof(UINT32)*KMER_LENGTH*NUMBER_OF_DNA);
    memset(pMostProbableKmer,0,sizeof(UINT32)*NUMBER_OF_DNA*2);
    
    //Random select kmers in DNA strings.
    for(i=0;i<NUMBER_OF_DNA;i++){
        srand(time(NULL));
        j = rand()%(SIZE_OF_DNA-KMER_LENGTH);
        memcpy(PROFILE_PATTERN++i*SIZE_OF_DNA,DNA_SRC+i*SIZE_OF_DNA+j,KMER_LENGTH);
    }
    
    best_score = score_mmotifs(PROFILE_PATTERN,KMER_LENGTH,NUMBER_OF_DNA);
    
    //Generate profile
    for(i=0;i<KMER_LENGTH;i++){
        a_cnt=c_cnt=g_cnt=t_cnt = 0;
        for(j=0;j<NUMBER_OF_DNA;j++){
            if(PROFILE_PATTERN[j*SIZE_OF_DNA+i] == 'A') a_cnt++;
            if(PROFILE_PATTERN[j*SIZE_OF_DNA+i] == 'C') c_cnt++;
            if(PROFILE_PATTERN[j*SIZE_OF_DNA+i] == 'G') g_cnt++;
            if(PROFILE_PATTERN[j*SIZE_OF_DNA+i] == 'T') t_cnt++;
        }
        PROFILE_PATTERN[i]               = a_cnt;
        PROFILE_PATTERN[SIZE_OF_DNA+i]   = c_cnt;
        PROFILE_PATTERN[2*SIZE_OF_DNA+i] = g_cnt;
        PROFILE_PATTERN[3*SIZE_OF_DNA+i] = t_cnt;
    }
        
    RandomizedMotifSearch(,,,,);
    
}
