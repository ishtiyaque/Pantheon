#pragma GCC diagnostic ignored "-Wnarrowing"

#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>
#include <sstream>
#include <cmath>
#include <ctime>
#include <stack>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include<cassert>

#include <iostream>
#include "seal/seal.h"

#include<globals.h>
#include "utils.hpp"
#include <openssl/sha.h>

using namespace std;
using namespace seal;

#define NUM_COL_THREAD (NUM_COL)
#define NUM_ROW_THREAD NUM_THREAD
#define NUM_PIR_THREAD 4

int TOTAL_MACHINE_THREAD = 6;

#define NUM_EXPANSION_THREAD (TOTAL_MACHINE_THREAD / NUM_COL)
#define NUM_EXPONENT_THREAD (TOTAL_MACHINE_THREAD / (NUM_COL_THREAD * NUM_ROW_THREAD))
// #define NUM_EXPONENT_THREAD (TOTAL_MACHINE_THREAD / (NUM_COL_THREAD ))

typedef  std::vector<seal::Ciphertext> PIRQuery;
typedef  seal::Ciphertext PIRResponse;

struct column_thread_arg {

    int col_id;
    int row_idx;
    Ciphertext *column_result;
    column_thread_arg(int c_id, int r_idx, Ciphertext *res) {
        col_id = c_id;
        row_idx = r_idx;
        column_result = res;
    }
};

struct mult_thread_arg {
    int id;
    int diff;
    Ciphertext *column_result;
    mult_thread_arg(int _id, int _diff, Ciphertext *_column_result) {
        id = _id;
        diff = _diff;
        column_result = _column_result;
    }
};

vector<vector<Plaintext>> db;

void sha256(const char *str, int len, unsigned char *dest);
void preallocate_memory();
void *expand_query(void *arg);
void *process_rows(void *arg);
void populate_db();
void *process_columns(void *arg);
void *multiply_columns(void *arg);
void *process_pir(void *arg);
void print_report();

void init_pir_params();
uint32_t get_next_power_of_two(uint32_t number);
uint32_t get_number_of_bits(uint64_t number);
void set_pir_db(std::vector<std::vector<uint64_t> > db);
void pir_encode_db(std::vector<std::vector<uint64_t>> db);
Ciphertext get_sum(Ciphertext *query, uint32_t start, uint32_t end);
vector<uint64_t> rotate_plain(std::vector<uint64_t> original, int index);

int NUM_ROW = 32;
int NUM_COL = 8;
int NUM_THREAD = 1;
vector<Plaintext> masks;
uint64_t number_of_items = 0;

SEALContext *context;

MemoryPoolHandle *column_pools;


Evaluator *evaluator;
BatchEncoder *batch_encoder;

Ciphertext *expanded_query;
Ciphertext *row_result;
Ciphertext *client_query_ct, *one_ct;
Ciphertext *server_query_ct;
Ciphertext *column_results;
Ciphertext *pir_results;

GaloisKeys galois_keys;
RelinKeys relin_keys;

seal::parms_id_type compact_pid;

pthread_t *query_expansion_thread;
int *expansion_thread_id;

pthread_t *row_process_thread;
int *row_thread_id;

pthread_t *pir_thread;
int *pir_thread_id;
//pthread_t *col_process_thread;
int *col_thread_id;

uint64_t *mult_time;
uint64_t *col_processing_time;
uint64_t *row_agg_time;


uint64_t query_gen_time, expansion_time, step1_time, step2_time, total_time;

/****PIR params ******************/

vector<Plaintext> pir_encoded_db;
uint32_t pir_num_obj;
uint32_t pir_obj_size;
uint32_t key_size; //  in bits
uint32_t pir_num_columns_per_obj;
uint32_t pir_num_query_ciphertext;
uint32_t pir_plain_bit_count;
uint32_t pir_db_rows;
/***********************************/

pthread_mutex_t print_lock = PTHREAD_MUTEX_INITIALIZER;

unsigned long *row_thread_processing_time_arr;

int main(int argc, char **argv)
{

    int option;
    const char *optstring = "n:k:s:t:";
    while ((option = getopt(argc, argv, optstring)) != -1)
    {
        switch (option)
        {
        case 'n':
            number_of_items = stoi(optarg);
            break;
        case 's':
            pir_obj_size = stoi(optarg);
            break;
        case 'k':
            key_size = stoi(optarg);
            break;
        case 't':
            NUM_THREAD = stoi(optarg);
            break;

        case '?':
            cout << "error optopt: " << optopt << endl;
            cout << "error opterr: " << opterr << endl;
            return 1;
        }
    }
    if(!number_of_items) {cout<<"Missing -n\n";return 0;}
    if(!NUM_THREAD) {cout<<"Missing -t\n";return 0;}
    if(!key_size) {cout<<"Missing -k\n"; return 0;}
    if(!pir_obj_size) {cout<<"Missing -s\n";return 0;} 

    NUM_COL = (int) ceil(key_size / (2.0 * PLAIN_BIT));
    NUM_ROW = (int) ceil(number_of_items / ((double)(N / 2)));
    init_pir_params();

    chrono::high_resolution_clock::time_point time_start, time_end, total_start, total_end;
    clock_t total_cpu_start, total_cpu_end, cpu_start, cpu_end;
    std::stringstream ss, qss;
    srand(time(NULL));

    mult_time = new uint64_t[NUM_ROW];
    col_processing_time = new uint64_t[NUM_ROW];
    row_agg_time = new uint64_t[NUM_ROW];

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(N);
    parms.set_coeff_modulus(CoeffModulus::Create(N, CT_PRIMES));


    parms.set_plain_modulus(PLAIN_MODULUS);

    context = new SEALContext(parms);
    auto pid = context->first_parms_id();
    uint64_t plain_modulus = parms.plain_modulus().value();
    KeyGenerator keygen(*context);
    SecretKey secret_key = keygen.secret_key();
    keygen.create_relin_keys(relin_keys);


    set<int> rotation_steps;  
    rotation_steps.insert(0);

    for (int i = N / (2 * NUM_COL); i < N / 2; i *= 2) {
        rotation_steps.insert(i);
    }

    for (int i = 1; i < (pir_num_columns_per_obj / 2); i *= 2)
    {
        rotation_steps.insert(-i);
    }
    keygen.create_galois_keys(vector<int>(rotation_steps.begin(), rotation_steps.end()), galois_keys);

    Encryptor encryptor(*context, secret_key);
    evaluator = new Evaluator(*context);
    Decryptor decryptor(*context, secret_key);

    column_pools = new MemoryPoolHandle[NUM_COL];
    for(int i = 0; i < NUM_COL; i++) {
        column_pools[i] = MemoryPoolHandle::New();
    }
    
    batch_encoder = new BatchEncoder(*context);
    size_t slot_count = batch_encoder->slot_count();
    size_t row_size = slot_count / 2;

    client_query_ct = new Ciphertext();
    server_query_ct = new Ciphertext();
    one_ct = new Ciphertext();

    vector<uint64_t> temp_mat;
    Plaintext temp_pt;

    vector<uint64_t> client_query_mat(N, 0ULL), result_mat, one_mat;
    for (int i = 0; i < N; i++)
    {
        one_mat.push_back(1);
    }
    Plaintext one_pt;
    batch_encoder->encode(one_mat, one_pt);

    encryptor.encrypt_symmetric(one_pt, *one_ct);

    for (int k = 0; k < MOD_SWITCH_COUNT; k++)
    {
        evaluator->mod_switch_to_next_inplace(*one_ct);
    }

    compact_pid = one_ct->parms_id();
    populate_db();
    vector<vector<uint64_t>> pir_db; //64 bit placeholder for 16 bit plaintext coefficients

    for(int i = 0; i < pir_num_obj; i++) {
        vector<uint64_t> v;
        for(int j = 0; j < (pir_obj_size / 2); j++) {  // 2 bytes each plaintxt slot
            v.push_back(rand() % PLAIN_MODULUS);
        }
        pir_db.push_back(v);
    }

    set_pir_db(pir_db);
    //cout << "DB population complete!" << endl;



    int desired_index = 100;
    int val = desired_index + 1;
    const char str[] = {val & 0xFF, (val >> 8) & 0xFF, (val >> 16) & 0xFF, (val >> 24) & 0xFF, 0};

	unsigned char hash[SHA256_DIGEST_LENGTH];
    sha256(str,4, hash);

    for(int i = 0; i < NUM_COL; i++) {
        for(int j = i * (N/(NUM_COL*2)); j < ((i + 1) * (N/(NUM_COL*2))); j++) {
            client_query_mat[j] = (uint64_t(hash[4*i]) << 8) + hash[4*i + 1];;
            client_query_mat[j + (N/2)] = (uint64_t(hash[4*i + 2]) << 8) + hash[4*i + 3];;
        }
    }

    row_thread_processing_time_arr = new unsigned long[NUM_ROW_THREAD];
    expanded_query = new Ciphertext[NUM_COL];
    pir_results = new Ciphertext[NUM_PIR_THREAD];

    row_result = new Ciphertext[NUM_ROW];

    row_thread_id = new int[NUM_ROW_THREAD];

    for (int i = 0; i < NUM_ROW_THREAD; i++)
    {
        row_thread_id[i] = i;
    }

    col_thread_id = new int[NUM_COL_THREAD];

    for (int i = 0; i < NUM_COL_THREAD; i++)
    {
        col_thread_id[i] = i;
    }

    Plaintext client_query_pt, result_pt;

    time_start = chrono::high_resolution_clock::now();

    batch_encoder->encode(client_query_mat, client_query_pt);
    Serializable<Ciphertext> ser_query = encryptor.encrypt_symmetric(client_query_pt);
    ser_query.save(qss);

    time_end = chrono::high_resolution_clock::now();
    query_gen_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();

    printf("query size (Byte): %lu\n", qss.str().size());


    query_expansion_thread = new pthread_t[NUM_COL];
    expansion_thread_id = new int[NUM_COL];

    for (int i = 0; i < NUM_COL; i++)
    {
        expansion_thread_id[i] = i;
    }

    row_process_thread = new pthread_t[NUM_ROW_THREAD];

    pir_thread = new pthread_t[NUM_PIR_THREAD];

    pir_thread_id = new int[NUM_PIR_THREAD];

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        pir_thread_id[i] = i;
    }

    for (int i = 0; i < NUM_COL; i++)
    {
        vector<uint64_t> mat(N, 0ULL);
        Plaintext pt;
        for (int j = i * (N / (2 * NUM_COL)); j < (i + 1) * (N / (2 * NUM_COL)); j++)
        {
            mat[j] = mat[j + (N / 2)] = 1;
        }
        batch_encoder->encode(mat, pt);
        evaluator->transform_to_ntt_inplace(pt, pid);
        masks.push_back(pt);
    }
    total_start = chrono::high_resolution_clock::now();
    total_cpu_start = clock();

    server_query_ct->load(*context, qss);

    total_start = chrono::high_resolution_clock::now();

    time_start = chrono::high_resolution_clock::now();
    cpu_start = clock();
    my_transform_to_ntt_inplace(*context, *server_query_ct, TOTAL_MACHINE_THREAD);
    for (int i = 0; i < NUM_COL; i++)
    {
        if (pthread_create(&(query_expansion_thread[i]), NULL, expand_query, (void *)&(expansion_thread_id[i])))
        {
            printf("Error creating expansion thread");
        }
    }

    for (int i = 0; i < NUM_COL; i++)
    {
        pthread_join(query_expansion_thread[i], NULL);
    }
    cpu_end = clock();
    time_end = chrono::high_resolution_clock::now();
    expansion_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
    auto query_expansion_cpu_time = double(cpu_end - cpu_start) / double(CLOCKS_PER_SEC);

    time_start = chrono::high_resolution_clock::now();
    cpu_start = clock();

    if(NUM_ROW_THREAD == 1) {
        process_rows((void *)&(row_thread_id[0]));
    }else {
        for (int i = 0; i < NUM_ROW_THREAD; i++)
        {
            if (pthread_create(&(row_process_thread[i]), NULL, process_rows, (void *)&(row_thread_id[i])))
            {
                printf("Error creating processing thread");
            }

        }

        for (int i = 0; i < NUM_ROW_THREAD; i++)
        {
            pthread_join(row_process_thread[i], NULL);
        }
    }

    time_end = chrono::high_resolution_clock::now();
    step1_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
   
    time_start = chrono::high_resolution_clock::now();

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        if (pthread_create(&(pir_thread[i]), NULL, process_pir, (void *)&(pir_thread_id[i])))
        {
            printf("Error creating PIR processing thread");
        }
    }

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        pthread_join(pir_thread[i], NULL);
    }
    for(int i = 1; i < NUM_PIR_THREAD; i++) {
        my_add_inplace(*context, pir_results[0], pir_results[i]);
    }
    cpu_end = clock();
    time_end = chrono::high_resolution_clock::now();
    step2_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
    auto processing_cpu_time = double(cpu_end - cpu_start) / double(CLOCKS_PER_SEC);
    total_time = (chrono::duration_cast<chrono::microseconds>(time_end - total_start)).count();

    Ciphertext final_result = pir_results[0];
    final_result.save(ss);

    total_cpu_end = clock();
    total_end = chrono::high_resolution_clock::now();

    auto latency = (chrono::duration_cast<chrono::microseconds>(total_end - total_start)).count();

    cout << "Result noise budget " << decryptor.invariant_noise_budget(final_result) << endl;
   
    decryptor.decrypt(final_result, result_pt);
    batch_encoder->decode(result_pt, result_mat);

    vector<uint64_t> decoded_response;
    decoded_response = rotate_plain(result_mat, desired_index % row_size);
    bool incorrect_result = false;
    for (int i = 0; i < pir_obj_size / 4 ; i++)
    {
        if ((pir_db[desired_index][i] != decoded_response[i]) || (pir_db[desired_index][i + pir_obj_size / 4] != decoded_response[i + N/2]))
        {
            incorrect_result = true;
            break;
        }
    }

    if (incorrect_result)
    {
        std::cout << "Result is incorrect!" << std::endl<<std::endl;
    }
    else
    {
        std::cout << "Result is correct!" << std::endl<<std::endl;
        print_report();
    }
    return 0;
}

void print_report() {
    
    cout<<"Query expansion time (ms): "<<expansion_time / 1000<<endl;
    cout<<"Equality check time (ms): "<<step1_time / 1000<<endl;
    cout<<"PIR time (ms): "<<step2_time / 1000<<endl;
    cout<<"Total processing time (ms): "<<total_time / 1000<<endl;
}

void populate_db()
{
    vector<vector<uint64_t> > mat_db;
    for(int i = 0;  i < NUM_ROW * NUM_COL; i++) {
        vector<uint64_t> v(N, 0ULL);
        mat_db.push_back(v);
    }
    unsigned char hash[SHA256_DIGEST_LENGTH];

    for (uint32_t row = 0; row < NUM_ROW * (N/2); row++)
    {
        uint32_t row_in_vector = row % (N/2);
        uint32_t val = row + 1;
        const char str[] = {val & 0xFF, (val >> 8) & 0xFF, (val >> 16) & 0xFF, (val >> 24) & 0xFF, 0};
        sha256(str,4, hash);
        for(int col = 0; col < NUM_COL; col++) {
            int vector_idx = (row / (N/2)) * NUM_COL + col;
            mat_db[vector_idx][row_in_vector] = (uint64_t(hash[4*col]) << 8) + hash[4*col + 1];
            mat_db[vector_idx][row_in_vector + (N/2)] = (uint64_t(hash[4*col + 2]) << 8) + hash[4*col + 3];

        }
    }

    for(int i = 0; i < NUM_ROW; i++) {
        vector<Plaintext> row_partition;
        for(int j = 0; j < NUM_COL; j++) {
            Plaintext pt;
            batch_encoder->encode(mat_db[i * NUM_COL + j], pt);
            row_partition.push_back(pt);
        }
        db.push_back(row_partition);
    }
    return;
}

void *expand_query(void *arg)
{
    int id = *((int *)arg);

    expanded_query[id] = *server_query_ct;
    my_multiply_plain_ntt(*context, expanded_query[id], masks[id], NUM_EXPANSION_THREAD);
    my_transform_from_ntt_inplace(*context, expanded_query[id], NUM_EXPANSION_THREAD);
    Ciphertext temp_ct;

    for (int i = N / (2 * NUM_COL); i < N / 2; i *= 2)
    {
        temp_ct = expanded_query[id];
        my_rotate_internal(*context, temp_ct, i, galois_keys, column_pools[id], NUM_EXPANSION_THREAD);
        my_add_inplace(*context, expanded_query[id], temp_ct);
    }
    return NULL;
}

void *process_rows(void *arg)
{
    int id = *((int *)arg);
    Ciphertext column_results[NUM_COL];
    Ciphertext pir_results[NUM_PIR_THREAD];
    vector<column_thread_arg> column_args;
    vector<mult_thread_arg> mult_args;

    for(int i = 0; i < NUM_COL_THREAD; i++) {
        column_args.push_back(column_thread_arg(i, id, column_results));
    }
    for(int i = 0; i < NUM_COL; i++) {
        mult_args.push_back(mult_thread_arg(i, 1, column_results));
    }
    int num_row_per_thread = NUM_ROW / NUM_ROW_THREAD;
    int start_idx = num_row_per_thread * id;
    int end_idx = start_idx + num_row_per_thread;

    pthread_t col_process_thread[NUM_COL_THREAD];
    pthread_t col_mult_thread[NUM_COL_THREAD];


    for(int row_idx = start_idx; row_idx < end_idx; row_idx++) {
        // time_start = chrono::high_resolution_clock::now();

        for (int i = 0; i < NUM_COL_THREAD; i++)
        {
            column_args[i].row_idx = row_idx;
            if (pthread_create(&(col_process_thread[i]), NULL, process_columns, (void *)&(column_args[i])))
            {
                printf("Error creating column processing thread");
            }
        }

        for (int i = 0; i < NUM_COL_THREAD; i++)
        {
            pthread_join(col_process_thread[i], NULL);
        }

        for(int diff = 2; diff <= NUM_COL; diff*=2) {

            for(int i = 0; i < mult_args.size(); i++) {
                mult_args[i].diff = diff;
            }
            for (int i = 0; i < NUM_COL; i+=diff)
            {
                if (pthread_create(&(col_mult_thread[i]), NULL, multiply_columns, (void *)&(mult_args[i])))
                {
                    printf("Error creating column processing thread");
                }
            }

            for (int i = 0; i < NUM_COL; i+=diff)
            {
                pthread_join(col_mult_thread[i], NULL);
            }

        }

        Ciphertext temp_ct = column_results[0];
        my_conjugate_internal(*context, temp_ct, galois_keys, column_pools[0], TOTAL_MACHINE_THREAD / NUM_ROW_THREAD);

        my_bfv_multiply(*context, column_results[0], temp_ct, column_pools[0], TOTAL_MACHINE_THREAD / NUM_ROW_THREAD);
        my_relinearize_internal(*context, column_results[0] ,relin_keys, 2, MemoryManager::GetPool(), TOTAL_MACHINE_THREAD / NUM_ROW_THREAD);
        my_transform_to_ntt_inplace(*context, column_results[0], TOTAL_MACHINE_THREAD );
        row_result[row_idx] = column_results[0];

    }

}

void *process_pir(void *arg) {
    int my_id = *((int *)arg);
    int column_per_thread = (pir_num_columns_per_obj / 2) / NUM_PIR_THREAD;
    int start_idx = my_id * column_per_thread;
    int end_idx = start_idx + column_per_thread - 1;
    pir_results[my_id] = get_sum(row_result, start_idx, end_idx);

    int mask = 1;
    while(mask <= start_idx) {
        if(start_idx & mask) {
            my_rotate_internal(*context, pir_results[my_id] , -mask, galois_keys, MemoryManager::GetPool(), TOTAL_MACHINE_THREAD/NUM_PIR_THREAD);
        }
        mask <<= 1;

    }
}


void *multiply_columns(void *arg) {
    mult_thread_arg mult_arg = *((mult_thread_arg *)arg);
    Ciphertext *column_results = mult_arg.column_result;
    int id = mult_arg.id;
    int diff = mult_arg.diff;
    int num_threads = TOTAL_MACHINE_THREAD / (NUM_COL / diff);

    my_bfv_multiply(*context, column_results[id], column_results[id + (diff/2)], column_pools[id], num_threads);
    my_relinearize_internal(*context, column_results[id] ,relin_keys, 2, column_pools[id], num_threads); 
}


void *process_columns(void *arg) {

    column_thread_arg col_arg = *((column_thread_arg *)arg);
    vector<unsigned long> exp_time;
    unsigned long exponent_time = 0;
    int num_col_per_thread = NUM_COL / NUM_COL_THREAD;
    int start_idx = num_col_per_thread * col_arg.col_id;
    int end_idx = start_idx + num_col_per_thread;
    for(int i = start_idx; i < end_idx; i++) {
        Ciphertext sub;
        Ciphertext prod;
        evaluator->sub_plain(expanded_query[i], db[col_arg.row_idx][i], sub);

        for (int k = 0; k < 16; k++){
            my_bfv_square(*context, sub, column_pools[i], NUM_EXPONENT_THREAD);
            my_relinearize_internal(*context, sub, relin_keys, 2, column_pools[i], NUM_EXPONENT_THREAD);
        }
        for(int k = 0; k < MOD_SWITCH_COUNT; k++) {
            my_mod_switch_scale_to_next(*context, sub, sub, column_pools[i], NUM_EXPONENT_THREAD);
        }
        evaluator->sub(*one_ct, sub, (col_arg.column_result)[i]);

    }
}



/*****************************PIR functions ************************************/
/******************************************************************************/

Ciphertext get_sum(Ciphertext *query, uint32_t start, uint32_t end)
{
    seal::Ciphertext result;

    if (start != end)
    {
        int count = (end - start) + 1;
        int next_power_of_two = get_next_power_of_two(count);
        int mid = next_power_of_two / 2;
        seal::Ciphertext left_sum = get_sum(query, start, start + mid - 1);
        seal::Ciphertext right_sum = get_sum(query, start + mid, end);
        my_rotate_internal(*context, right_sum, -mid, galois_keys, column_pools[0], TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
        my_add_inplace(*context, left_sum, right_sum);
        return left_sum;
    }
    else
    {

        seal::Ciphertext column_sum = query[0];
        seal::Ciphertext temp_ct;
        my_multiply_plain_ntt(*context, column_sum, pir_encoded_db[pir_num_query_ciphertext * start], TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);

        for (int j = 1; j < pir_num_query_ciphertext; j++)
        {
            temp_ct = query[j];
            my_multiply_plain_ntt(*context, temp_ct, pir_encoded_db[pir_num_query_ciphertext * start + j], TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
            my_add_inplace(*context, column_sum, temp_ct);
        }
        my_transform_from_ntt_inplace(*context, column_sum, TOTAL_MACHINE_THREAD / NUM_PIR_THREAD);
        return column_sum;
    }
}

uint32_t get_next_power_of_two(uint32_t number)
{
    if (!(number & (number - 1)))
    {
        return number;
    }

    uint32_t number_of_bits = get_number_of_bits(number);
    return (1 << number_of_bits);
}


uint32_t get_number_of_bits(uint64_t number)
{
    uint32_t count = 0;
    while (number)
    {
        count++;
        number /= 2;
    }
    return count;
}


void set_pir_db(std::vector<std::vector<uint64_t> > db)
{
    assert(db.size() == pir_num_obj);
    std::vector<std::vector<uint64_t> > extended_db(pir_db_rows);
    for(int i = 0; i < pir_db_rows; i++) {
        extended_db[i] = std::vector<uint64_t>(N, 1ULL);
    }
    int row_size = N/2;

    for(int i = 0; i < pir_num_obj;i++) {
        std::vector<uint64_t> temp = db[i];

        int row = (i / row_size);
        int col = (i % row_size);
        for (int j = 0; j < pir_num_columns_per_obj / 2; j++)
        {
            extended_db[row][col] = temp[j];
            extended_db[row][col+row_size] = temp[j+(pir_num_columns_per_obj / 2)];
            row += pir_num_query_ciphertext;
        }

    }   
    pir_encode_db(extended_db);
    return;
}

void pir_encode_db(std::vector<std::vector<uint64_t>> db)
{
    pir_encoded_db = std::vector<seal::Plaintext>(db.size());
    for (int i = 0; i < db.size(); i++)
    {
        batch_encoder->encode(db[i], pir_encoded_db[i]);
        evaluator->transform_to_ntt_inplace(pir_encoded_db[i], compact_pid);

    }
}

vector<uint64_t> rotate_plain(std::vector<uint64_t> original, int index)
{
    int sz = original.size();
    int row_count = sz / 2;
    std::vector<uint64_t> result(sz);
    for (int i = 0; i < row_count; i++)
    {
        result[i] = original[(index + i) % row_count];
        result[row_count + i] = original[row_count + ((index + i) % row_count)];
    }

    return result;
}

void init_pir_params()
{
    pir_plain_bit_count = 16;
    pir_num_obj = ((N/2) * NUM_ROW);
    pir_num_query_ciphertext = ceil(pir_num_obj / (double)(N/2));
    pir_num_columns_per_obj = 2 * (ceil(((pir_obj_size/2) * 8) / (float)(pir_plain_bit_count)));
    pir_db_rows = ceil(pir_num_obj / (double)N) * pir_num_columns_per_obj;

    return;

}

void sha256(const char *str, int len, unsigned char *dest)
{    
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, str, len);
    SHA256_Final(dest, &sha256);
}