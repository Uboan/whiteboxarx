#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdbool.h>
#include "white_box_backend.c"
# define LLONG_MAX __LONG_LONG_MAX__

/* Maximum value an `unsigned long long int' can hold.  (Minimum is 0).  */
# undef ULLONG_MAX
# define ULLONG_MAX (LLONG_MAX * 2ULL + 1ULL)
#define rotate_left(x,n) (((x) >> (WORD_SIZE - (n))) | ((x) << (n)))
#define rotate_right(x,n) (((x) << (WORD_SIZE - (n))) | ((x) >> (n)))



uint64_t mod_mask = ULLONG_MAX >> (64 - WORD_SIZE);

WORD_TYPE mod_substract(WORD_TYPE a,WORD_TYPE b){
    uint64_t tmp,tmp2;
    tmp = a & mod_mask;
    tmp2 = b & mod_mask;
   // printf("\ntmp = %lld\n",tmp);
   // printf("mod_mask = %lld\n",mod_mask);
   // printf("tmp2 = %lld\n",tmp2);
    tmp2 = (mod_mask+1)-tmp2;
   // printf("inversed tmp2 = %lld\n",tmp2);
    tmp2 = tmp+tmp2;
   // printf("tmp + tmp2 = %lld\n",tmp2);
    if(tmp2>mod_mask){
        tmp2=tmp2 - mod_mask-1;
    }
    a =  tmp2;
   // printf("a=%lld", a);
    return a;
}

bool Recover_KEY(WORD_TYPE *guess_key,WORD_TYPE *recovered_key){
    if(guess_key==NULL&&recovered_key==NULL)return false;
    
    recovered_key[3]=guess_key[0];//k0
    WORD_TYPE tmp;
    tmp = guess_key[1]^0^rotate_left(guess_key[0],2);
    recovered_key[2] = mod_substract(tmp,guess_key[0]);
    recovered_key[2]=rotate_left(recovered_key[2],7);


    tmp = guess_key[2]^1^rotate_left(guess_key[1],2);
    recovered_key[1] = mod_substract(tmp,guess_key[1]);
    recovered_key[1]=rotate_left(recovered_key[1],7);

    tmp = guess_key[3]^2^rotate_left(guess_key[2],2);
    recovered_key[0] = mod_substract(tmp,guess_key[2]);
    recovered_key[0]=rotate_left(recovered_key[0],7);
    return true;
} 
#define ROTATE_RIGHT(x, pos) ((x >> pos) | (x << (WORD_SIZE - pos)))
#define ROTATE_LEFT(x, pos) ((x << pos) | (x >> (WORD_SIZE - pos)))
#define KEY_WORDS  4
#define ROUND(k, x, y) ( \
    x = ROTATE_RIGHT(x, 7), \
    x += y, \
    x ^= k, \
    y = ROTATE_LEFT(y, 2), \
    y ^= x \
)
void key_expansion(WORD_TYPE key[KEY_WORDS], WORD_TYPE k[ROUNDS]) {
    k[0] = key[KEY_WORDS - 1];

    WORD_TYPE l[KEY_WORDS - 1 + ROUNDS - 1];
    for (size_t i = 1; i < KEY_WORDS; i++) {
        l[i - 1] = key[KEY_WORDS - 1 - i];
    }

    for (size_t i = 0; i < ROUNDS - 1; i++) {
        l[KEY_WORDS - 2 + i + 1] = l[i];
        k[i + 1] = k[i];
        ROUND(i, l[KEY_WORDS - 2 + i + 1], k[i + 1]);
    }
}



int main(int argc, char *argv[]){
    WORD_TYPE GK[4],KEY[4],k[ROUNDS];
    if(argc<5){
        return -1;
    }else{
        sscanf(argv[1],"%" WORD_IN_TYPE,&GK[0]);
        sscanf(argv[2],"%" WORD_IN_TYPE,&GK[1]);
        sscanf(argv[3],"%" WORD_IN_TYPE,&GK[2]);
        sscanf(argv[4],"%" WORD_IN_TYPE,&GK[3]);
        printf("%" WORD_OUT_TYPE " %" WORD_OUT_TYPE " %" WORD_OUT_TYPE " %" WORD_OUT_TYPE "\n", GK[0],GK[1],GK[2],GK[3]);
        Recover_KEY(GK,KEY);
       
        printf("\n%d\n",WORD_SIZE);
        printf("%" WORD_OUT_TYPE " %" WORD_OUT_TYPE " %" WORD_OUT_TYPE " %" WORD_OUT_TYPE "\n", KEY[0], KEY[1], KEY[2], KEY[3]);
   
        // key_expansion(GK,k);
        // for(int i=0;i<ROUNDS;i++){
        //     printf("%" WORD_OUT_TYPE " ",k[i]);
            
        // }
        
 
        return 1;

    }



}