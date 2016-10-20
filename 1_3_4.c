/*
 * Gibbs Sampling
 * Not the Profile-randomly generated k-mer 
 * Delete the latest "\n" from input file
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#define DEBUG

typedef unsigned char    UINT8;
typedef unsigned int     UINT32;


#define MAX_ARRAY_SIZE   1000000
#define MAX_PROCESS_TIME 60

#define SCAN_KMER         0
#define SCAN_NUMOFDNA     1
#define SCAN_DNA          2
#define SCAN_PROCESS_TIME 3


typedef struct {
    UINT8 *src;
    UINT32 size_of_each_dna;
    UINT32 kmer;
    UINT32 num_of_dna;
    UINT32 total_size;
}DNA_INFO_t;

UINT32 alpha[256] = {0};

/*
 * Read and parse data.
 */
bool read_test_data(const char* file_name,DNA_INFO_t *dna_info,UINT32 *N,int maxsize) 
{    
    UINT32 phase = SCAN_KMER,count=0;
    UINT8 *src = dna_info->src;
    UINT8 c;
    FILE* file = fopen (file_name, "r");
    
    if (!file) return false;
    if (!src)  return false;
    
    *N = dna_info->total_size = dna_info->num_of_dna = dna_info->kmer = 0;
    
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
            if(c == ' '){
                phase = SCAN_PROCESS_TIME;
                continue;
            }

            dna_info->num_of_dna *= 10;
            dna_info->num_of_dna += c-'0';
        }
        
        if(SCAN_PROCESS_TIME == phase){
            c = fgetc(file);
            if((c == 0xd) || (c == 0xa)){
                if(0xa == c)
                    phase = SCAN_DNA;
                continue;
            }
            
            *N *= 10;
            *N += c-'0';
        }

    } while (!feof (file) && count < maxsize);
    
    dna_info->size_of_each_dna = dna_info->total_size/dna_info->num_of_dna;
    
    fclose (file);
    
    return true;
}



/*
 * dna_info: 
 * profile_pattern:
 */
double caculate_pr(UINT8 *dna,UINT32 *profile_pattern,UINT32 kmer)
{
    UINT32 i;
    double multiple = 1;
    double j;
    
    for(i=0,multiple=1;i<kmer;i++){
        
        j = profile_pattern[i]+profile_pattern[i+kmer]+profile_pattern[i+kmer*2]+profile_pattern[i+kmer*3];
        j = profile_pattern[alpha[dna[i]]*kmer+i]/j;
                    
        multiple*=j;
    }

    return multiple;
}


/*
 * Score the matrix motifs.
 * score_dna
 * dna_info
 */
UINT32 score_mmotifs(UINT8 *score_dna,DNA_INFO_t *dna_info)
{
    UINT32 i,j,l;
    UINT32 sum_A,sum_C,sum_G,sum_T,max_NU=0,sum_column=0,sum=0;
    UINT32 kmer = dna_info->kmer;
    UINT32 num_dna = dna_info->num_of_dna;
    UINT32 size_dna = dna_info->size_of_each_dna;

    for(i=0;i<kmer;i++){        
        for(sum_A=sum_C=sum_G=sum_T=j=0;j<num_dna;j++){           
            if(score_dna[j*kmer+i] == 'A') sum_A++;
            if(score_dna[j*kmer+i] == 'C') sum_C++;
            if(score_dna[j*kmer+i] == 'G') sum_G++;
            if(score_dna[j*kmer+i] == 'T') sum_T++;
        }
        sum_column = sum_A+sum_C+sum_G+sum_T;
        l = (sum_A>sum_C)?sum_A:sum_C;
        l = (sum_G>l)?sum_G:l;
        l = (sum_T>l)?sum_T:l;
        sum_column -= l;
        sum += sum_column;
    }

    return sum;
}
/*
 * Generate motifs profile.
 * motifs:
 * profile_pattern:
 * dna_info:
 * ith: except ith sequence.
 */
void gen_profile(UINT8 *motifs,UINT32 *profile_pattern,DNA_INFO_t *dna_info,UINT32 ith_except)
{
    UINT32 i,j;
    UINT32 a_cnt,c_cnt,g_cnt,t_cnt;
    
    //Generate profile
    for(i=0;i<dna_info->kmer;i++){  
        for(a_cnt=c_cnt=g_cnt=t_cnt=j=0;j<dna_info->num_of_dna;j++){
            if(ith_except==j) continue;
            if(motifs[j*dna_info->kmer+i] == 'A') a_cnt++;
            if(motifs[j*dna_info->kmer+i] == 'C') c_cnt++;
            if(motifs[j*dna_info->kmer+i] == 'G') g_cnt++;
            if(motifs[j*dna_info->kmer+i] == 'T') t_cnt++;
        }
        profile_pattern[i]                  = a_cnt+1;
        profile_pattern[dna_info->kmer+i]   = c_cnt+1;
        profile_pattern[2*dna_info->kmer+i] = g_cnt+1;
        profile_pattern[3*dna_info->kmer+i] = t_cnt+1;
    }
}


/*
 * Limit to 1000 times
 * dna_info: .
 * best_motifs: The first randon select motifs.
 * N: The number of the random process.
 */
UINT32 GibbsSampler(DNA_INFO_t *dna_info,UINT32 *profile_pattern,UINT8 *best_motifs,UINT32 N) 
{
    UINT32 score=0;
    UINT32 i,j,ith_except;
    UINT8  *dna_src = dna_info->src;
    const UINT32 kmer = dna_info->kmer;
    const UINT32 num_dna = dna_info->num_of_dna;
    const UINT32 size_dna = dna_info->size_of_each_dna;
    UINT32 best_score;
    UINT8  *motifs = malloc(kmer*num_dna);    
    double m_prob,c_prob;

    //Random select motifs from DNA strings.
    for(i=0;i<num_dna;i++){
        j = rand()%(size_dna-kmer);
        memcpy(best_motifs+i*kmer,dna_src+i*size_dna+j,kmer);
    }
    
    memcpy(motifs,best_motifs,kmer*num_dna);
    
    // Get the score from first random motifs.
    best_score = score_mmotifs(best_motifs,dna_info);
    
    for(i=0;i<N;i++){
        // i ¡ö Random(t)
        ith_except = rand()%num_dna;
        // Profile ¡ö profile matrix constructed from all strings in Motifs except for Motifi
        gen_profile(motifs,profile_pattern,dna_info,ith_except);
        
        for(j=0,m_prob=0;j<=(size_dna-kmer);j++){              
                  c_prob = caculate_pr(dna_src+j+ith_except*size_dna,profile_pattern,kmer);
                  if(c_prob>m_prob){
                      m_prob = c_prob;
                      //Update the max probabity position
                      memcpy(motifs+ith_except*kmer,dna_src+j+ith_except*size_dna,kmer);
                  }
        }
        //Score
        score = score_mmotifs(motifs,dna_info);
        if(score<best_score){
            best_score = score;
            memcpy(best_motifs,motifs,kmer*num_dna);
        }
    }
    return best_score;
}

void main()
{
    UINT32 i=0,j=0,N=0;
    UINT32 pr_multi,max_pr=0;
    UINT32 *profile_pattern;
    DNA_INFO_t dna_info;
    UINT8  *best_motifs,*motifs;
    UINT32 score = 0xffffffff;
    time_t t;
        
    dna_info.src = malloc(MAX_ARRAY_SIZE);
    if(!read_test_data("dataset_158_9.txt",&dna_info,&N,MAX_ARRAY_SIZE))
        return;
#ifdef DEBUG
    printf("num_d:%d dna_total:%d dna_sz:%d kmer:%d process%d \n",dna_info.num_of_dna,dna_info.total_size,dna_info.size_of_each_dna,dna_info.kmer,N);
#endif
    // Init hash index.
    alpha['A'] = 0;
    alpha['C'] = 1;
    alpha['G'] = 2;
    alpha['T'] = 3;

    //Init profile table
    profile_pattern = malloc(sizeof(UINT32)*dna_info.kmer*4);
    memset(profile_pattern,0,sizeof(UINT32)*dna_info.kmer*4);    
    
    best_motifs = malloc(dna_info.kmer*dna_info.num_of_dna);
    motifs = malloc(dna_info.kmer*dna_info.num_of_dna);
    //Initialize random info.
    srand((unsigned) time(&t));
    
    for(i=0,score=0xffffffff;i<MAX_PROCESS_TIME;i++){
        j = GibbsSampler(&dna_info,profile_pattern,motifs,N);
        if(j<score){
            score = j;
            memcpy(best_motifs,motifs,dna_info.kmer*dna_info.num_of_dna);
        }
    }

#ifdef DEBUG
    printf("score %d \n",score);
#endif

    for(i=0;i<dna_info.num_of_dna;i++){
        for(j=0;j<dna_info.kmer;j++){
            printf("%c",best_motifs[i*dna_info.kmer+j]);
        }
        printf("\n"); 
    }
    
    free(best_motifs);
}
