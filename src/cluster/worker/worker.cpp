#include<globals.h>
#include<queue>
#include <semaphore.h>

using namespace std::chrono;
using namespace std;
using namespace seal;


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

int NUM_TOTAL_COL = 0;
int NUM_ROW = 0;
int NUM_COL_THIS_WORKER = 0;
int NUM_THREAD = 1;
vector<Plaintext> masks;

int WORKER_ID = -1;
int NUM_TOTAL_WORKER = 0;
int GROUP_SIZE = 0;
bool GROUP_LEADER = false;
int ID_IN_GROUP;
int NUM_GROUP;

#define NUM_COL_THREAD (NUM_COL_THIS_WORKER)
#define NUM_ROW_THREAD NUM_THREAD
#define NUM_PIR_THREAD 4

int TOTAL_MACHINE_THREAD = 24;

#define NUM_EXPANSION_THREAD (TOTAL_MACHINE_THREAD / NUM_COL_THIS_WORKER)
//#define NUM_EXPONENT_THREAD (TOTAL_MACHINE_THREAD / (NUM_COL_THREAD * NUM_ROW_THREAD))
#define NUM_EXPONENT_THREAD (TOTAL_MACHINE_THREAD / (NUM_COL_THREAD ))


/****PIR params ******************/

vector<Plaintext> pir_encoded_db;
uint32_t pir_num_obj;
uint32_t pir_obj_size;
uint32_t pir_num_columns_per_obj;
uint32_t pir_num_query_ciphertext;
uint32_t pir_plain_bit_count;
uint32_t pir_db_rows;
/***********************************/

SEALContext *context;

MemoryPoolHandle *column_pools;


Evaluator *evaluator;
BatchEncoder *batch_encoder;

Ciphertext *expanded_query;
Ciphertext *row_result;
Ciphertext *one_ct;
Ciphertext *server_query_ct;
Ciphertext *column_results;
Ciphertext *pir_results;

vector<vector<Plaintext>> db;
// vector<Plaintext> mult_pt;
queue<Ciphertext> group_response_queue;

vector<uint64_t> worker_response;

string *worker_ip;
string serialized_query;

GaloisKeys *galois_keys;
RelinKeys *relin_keys;
clock_t total_cpu_start_time, total_cpu_stop_time;

seal::parms_id_type compact_pid;

pthread_t *query_expansion_thread;
int *expansion_thread_id;

pthread_t *row_process_thread;  // Need to remove
int *row_thread_id;

//pthread_t *col_process_thread;
// int *col_thread_id;

pthread_t *pir_thread;
int *pir_thread_id;

uint64_t *mult_time;
uint64_t *col_processing_time;
uint64_t *row_agg_time;

pthread_cond_t sent_response_cond = PTHREAD_COND_INITIALIZER;
pthread_mutex_t sent_response_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t group_response_lock = PTHREAD_MUTEX_INITIALIZER;


sem_t group_response_sem;

uint64_t query_rcv_timestamp, response_send_timestamp;
uint64_t pir_start_timestamp, pir_end_timestamp;
uint64_t processing_start_timestamp, processing_end_timestamp;
uint64_t response_send_start_timestamp,response_send_end_timestamp;
uint64_t processing_clock_time;
double processing_cpu_time;

chrono::high_resolution_clock::time_point total_start, total_end;

void readWorkerIp();
vector<uint64_t> sendQuery(string _serialized_query);
void sendKeys(string _serialized_gal_key, string _serialized_relin_key, string _serialized_one_ct);
void printReport();
void preallocate_memory();
void *expand_query(void *arg);
void *process_rows(void *arg);
void populate_db();
void *process_columns(void *arg);
void *multiply_columns(void *arg);
void sendGroupResponse(vector<uint64_t> ct);
void init_pir_params();
void populate_pir_db();
void *process_pir(void *arg);
uint32_t get_next_power_of_two(uint32_t number);
uint32_t get_number_of_bits(uint64_t number);
Ciphertext get_sum(Ciphertext *query, uint32_t start, uint32_t end);
void set_pir_db(std::vector<std::vector<uint64_t> > db);
void pir_encode_db(std::vector<std::vector<uint64_t>> db);

rpc::client **group_client;
rpc::server *group_server;

int main(int argc, char *argv[]) {

    uint64_t number_of_items = 0;
    chrono::high_resolution_clock::time_point time_start, time_end;

    int option;
    const char *optstring = "n:t:i:w:c:g:s:";
    while ((option = getopt(argc, argv, optstring)) != -1) {
	switch (option) {

        case 'n':
            number_of_items = stoi(optarg);
            break;
        case 't':
            NUM_THREAD = stoi(optarg);
            break;
        case 'i':
            WORKER_ID = stoi(optarg);
            break;
        case 'w':
            NUM_TOTAL_WORKER = stoi(optarg);
            break; 
        case 'c':
            NUM_TOTAL_COL = stoi(optarg);
            break; 
        case 'g':
            GROUP_SIZE = stoi(optarg);
            break;  
        case 's':
            pir_obj_size = stoi(optarg);
            break;  
        case '?':
            cout<<"error optopt: "<<optopt<<endl;
            cout<<"error opterr: "<<opterr<<endl;
            return 1;
	    }
    }
    if(!number_of_items) {cout<<"Missing -n\n";return 0;}
    if(!NUM_THREAD) {cout<<"Missing -t\n";return 0;}
    if(!NUM_TOTAL_WORKER) {cout<<"Missing -w\n";return 0;}
    if(WORKER_ID < 0) {cout<<"Missing -i\n";return 0;}
    if(!GROUP_SIZE) {cout<<"Missing -g\n"; return 0;}
    if(!NUM_TOTAL_COL) {cout<<"Missing -c\n"; return 0;}
    if(!pir_obj_size) {cout<<"Missing -s\n";return 0;} 

    NUM_GROUP = (NUM_TOTAL_WORKER / GROUP_SIZE);

    NUM_ROW = (int) ceil(number_of_items / ((double)N / 2));
    NUM_ROW = NUM_ROW / NUM_GROUP;
    ID_IN_GROUP = WORKER_ID % GROUP_SIZE;
    if(ID_IN_GROUP == 0) {
        GROUP_LEADER = true;
    }
    NUM_COL_THIS_WORKER = NUM_TOTAL_COL / GROUP_SIZE;

    init_pir_params();
    readWorkerIp();

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(N);
    //parms.set_coeff_modulus(CoeffModulus::BFVDefault(N));
    parms.set_coeff_modulus(CoeffModulus::Create(N, CT_PRIMES));

    parms.set_plain_modulus(PLAIN_MODULUS);
    context = new SEALContext(parms);

    auto pid = context->first_parms_id();
    uint64_t plain_modulus = parms.plain_modulus().value();
    //print_parameters(*context);
    //cout << "Maximum size of co-efficient modulus: " << CoeffModulus::MaxBitCount(N) << endl;
    //auto qualifiers = context->first_context_data()->qualifiers();
    //cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;
    //cout << "Fast plain lift enabled: " << qualifiers.using_fast_plain_lift << endl;

    evaluator = new Evaluator(*context);

    column_pools = new MemoryPoolHandle[NUM_COL_THIS_WORKER];
    for(int i = 0; i < NUM_COL_THIS_WORKER; i++) {
        column_pools[i] = MemoryPoolHandle::New();
    }

    preallocate_memory();

    galois_keys = new GaloisKeys();
    relin_keys = new RelinKeys();
    one_ct = new Ciphertext();

    batch_encoder = new BatchEncoder(*context);
    size_t slot_count = batch_encoder->slot_count();
    size_t row_size = slot_count / 2;

    server_query_ct = new Ciphertext();

    populate_db();

    // Allocating memory
    for(int i = 0 ; i < SMALL_COEFF_COUNT; i++) {  
        worker_response.push_back(1ULL);
    }

    sem_init(&group_response_sem, 0, 0);

    for (int i = (ID_IN_GROUP * NUM_COL_THIS_WORKER); i < ((ID_IN_GROUP + 1) * NUM_COL_THIS_WORKER); i++)
    {
        vector<uint64_t> mat(N, 0ULL);
        Plaintext pt;
        for (int j = i * (N / (2 * NUM_TOTAL_COL)); j < (i + 1) * (N / (2 * NUM_TOTAL_COL)); j++)
        {
            mat[j] = mat[j + (N / 2)] = 1;
        }
        batch_encoder->encode(mat, pt);
        evaluator->transform_to_ntt_inplace(pt, pid);
        masks.push_back(pt);
    }

    query_expansion_thread = new pthread_t[NUM_COL_THIS_WORKER];
    expansion_thread_id = new int[NUM_COL_THIS_WORKER];

    for (int i = 0; i < NUM_COL_THIS_WORKER; i++)
    {
        expansion_thread_id[i] = i;
    }
    expanded_query = new Ciphertext[NUM_COL_THIS_WORKER];

    row_thread_id = new int[NUM_ROW_THREAD];

    for (int i = 0; i < NUM_ROW_THREAD; i++)
    {
        row_thread_id[i] = i;
    }

    pir_results = new Ciphertext[NUM_PIR_THREAD];
    pir_thread = new pthread_t[NUM_PIR_THREAD];

    pir_thread_id = new int[NUM_PIR_THREAD];

    for (int i = 0; i < NUM_PIR_THREAD; i++)
    {
        pir_thread_id[i] = i;
    }

    // col_thread_id = new int[NUM_COL_THREAD];

    // for (int i = 0; i < NUM_COL_THREAD; i++)
    // {
    //     col_thread_id[i] = i;
    // }

    row_result = new Ciphertext[NUM_ROW];

    rpc::server query_server(string(worker_ip[WORKER_ID]), WORKER_PORT + WORKER_ID);
    query_server.bind("sendKeys", sendKeys);
    query_server.bind("sendQuery", sendQuery);

    if(GROUP_LEADER) {
        query_server.bind("sendGroupResponse", sendGroupResponse);
        query_server.async_run(2);

    } else {
        query_server.async_run(1);
    }

    //cout<<"started query server"<<endl;   

    pthread_mutex_lock(&sent_response_lock);
    pthread_cond_wait(&sent_response_cond, &sent_response_lock);
    pthread_mutex_unlock(&sent_response_lock);

    total_end = chrono::high_resolution_clock::now();
    total_cpu_stop_time = clock();

    sleep(3);
    processing_clock_time  = chrono::duration_cast<chrono::microseconds>(total_end - total_start).count();
    processing_cpu_time = (((double)total_cpu_stop_time-total_cpu_start_time)/CLOCKS_PER_SEC);

    printReport();

    return 0;
}

vector<uint64_t> sendQuery(string _serialized_query) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    // cout<<"received query"<<endl;

    total_cpu_start_time = clock();
    query_rcv_timestamp = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());
    total_start = chrono::high_resolution_clock::now();

    if(GROUP_LEADER) {
        for(int i = 0; i < (GROUP_SIZE - 1); i++) {
            group_client[i]->async_call("sendQuery", _serialized_query);
        }
    }
    stringstream qss(_serialized_query);
    server_query_ct->load(*context, qss);

    my_transform_to_ntt_inplace(*context, *server_query_ct, TOTAL_MACHINE_THREAD);
    for (int i = 0; i < NUM_COL_THIS_WORKER; i++)
    {
        if (pthread_create(&(query_expansion_thread[i]), NULL, expand_query, (void *)&(expansion_thread_id[i])))
        {
            printf("Error creating expansion thread");
        }
    }

    for (int i = 0; i < NUM_COL_THIS_WORKER; i++)
    {
        pthread_join(query_expansion_thread[i], NULL);
    }
    //time_end = chrono::high_resolution_clock::now();
    //auto query_expansion_time = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
    //auto query_expansion_cpu_time = double(cpu_end - cpu_start) / double(CLOCKS_PER_SEC);

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
    if(!GROUP_LEADER) {
        pthread_mutex_lock(&sent_response_lock);
        pthread_cond_broadcast(&sent_response_cond);
        pthread_mutex_unlock(&sent_response_lock);
        
        response_send_timestamp = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());
        return vector<uint64_t>();
    }

    pir_start_timestamp = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());

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
        // evaluator->add_inplace(pir_results[0], pir_results[i]);
        my_add_inplace(*context, pir_results[0], pir_results[i]);
    }

    std::copy(pir_results[0].data(), pir_results[0].data() + SMALL_COEFF_COUNT, worker_response.begin());

    // PIRResponse pir_response = get_sum(row_result, 0, pir_num_columns_per_obj / 2 - 1);
    pir_end_timestamp = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());

    // std::copy(pir_response.data(), pir_response.data() + SMALL_COEFF_COUNT, worker_response.begin());
    pthread_mutex_lock(&sent_response_lock);
    pthread_cond_broadcast(&sent_response_cond);
    pthread_mutex_unlock(&sent_response_lock);
    
    response_send_timestamp = (chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count());
    return worker_response;
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
            // evaluator->rotate_rows_inplace(pir_results[my_id], -mask, *galois_keys);
            my_rotate_internal(*context, pir_results[my_id] , -mask, *galois_keys, MemoryManager::GetPool(), TOTAL_MACHINE_THREAD/NUM_PIR_THREAD);

        }
        mask <<= 1;

    }
}

void sendKeys(string _serialized_gal_key, string _serialized_relin_key, string _serialized_one_ct) {

    // cout<<"received keys"<<endl;
    stringstream gkss(_serialized_gal_key);
    stringstream rkss(_serialized_relin_key);
    stringstream oss(_serialized_one_ct);
    galois_keys->load(*context,gkss);

    relin_keys->load(*context, rkss);

    one_ct->load(*context, oss);

    for (int k = 0; k < MOD_SWITCH_COUNT; k++)
    {
        evaluator->mod_switch_to_next_inplace(*one_ct);
    }

    compact_pid = one_ct->parms_id();
    if(!GROUP_LEADER) {
        group_client = new rpc::client*[1];
        int leader_id = (WORKER_ID / GROUP_SIZE) * GROUP_SIZE;
        group_client[0] = new rpc::client(worker_ip[leader_id], WORKER_PORT + leader_id);
        return;
    }

    group_client = new rpc::client*[GROUP_SIZE - 1];
    for(int i = 0; i < (GROUP_SIZE - 1); i++) {
        int id = WORKER_ID + 1 + i;
        group_client[i] = new rpc::client(worker_ip[id], WORKER_PORT + id);
    }

    // for (int i = 0; i < NUM_ROW; i++)
    // {
    //     vector<uint64_t> mat(N, ((WORKER_ID / GROUP_SIZE) * NUM_ROW) + i + 1);
    //     Plaintext pt;
    //     batch_encoder->encode(mat, pt);
    //     evaluator->transform_to_ntt_inplace(pt,compact_pid);
    //     mult_pt.push_back(pt);
    // }
    if(GROUP_LEADER) {
        populate_pir_db();
    }

    return;

}


void readWorkerIp() {
    worker_ip = new string[NUM_TOTAL_WORKER];

    fstream fs;
	fs.open("../common/worker_ip.txt", ios::in);
	if (!fs) {
		cout << "Missing worker_ip.txt\n";
        exit(1);
	}

    for(int i = 0; i< NUM_TOTAL_WORKER;  i++) {
        fs>>worker_ip[i];
    }

    return;

}

void printReport() {
    cout<<"query receive timestamp"<<endl<<query_rcv_timestamp<<endl;

    cout<<"\nresponse send timestamp"<<endl<<response_send_timestamp<<endl;

    cout<<"\nTotal processing clock time"<<endl<<processing_clock_time<<endl;

    cout<<"\ntotal cpu time (sec) "<<endl<<processing_cpu_time<<endl;

    cout<<"\nPIR clock time "<<endl<<(pir_end_timestamp - pir_start_timestamp)<<endl;

    return;
}

void populate_db()
{
    for (int i = 0; i < NUM_ROW; i++)
    {
        vector<Plaintext> row;
        for (int j = 0; j < NUM_COL_THIS_WORKER; j++)
        {
            vector<uint64_t> mat;
            Plaintext pt;
            for (int k = 0; k < N; k++)
            {
                // mat.push_back((j * NUM_ROW + i * NUM_COL_THIS_WORKER + k) % PLAIN_MODULUS);
                mat.push_back(rand() % PLAIN_MODULUS);
            }
            batch_encoder->encode(mat, pt);
            row.push_back(pt);
        }
        db.push_back(row);
    }
    return;
}

void *expand_query(void *arg)
{
    int id = *((int *)arg);
    // evaluator->multiply_plain(*server_query_ct, masks[id], expanded_query[id]);
    // evaluator->transform_from_ntt_inplace(expanded_query[id]);

    expanded_query[id] = *server_query_ct;
    my_multiply_plain_ntt(*context, expanded_query[id], masks[id], NUM_EXPANSION_THREAD);
    my_transform_from_ntt_inplace(*context, expanded_query[id], NUM_EXPANSION_THREAD);
    Ciphertext temp_ct;

    for (int i = N / (2 * NUM_COL_THIS_WORKER); i < N / 2; i *= 2)
    {
        //evaluator->rotate_rows(expanded_query[id], i, *galois_keys, temp_ct);
        temp_ct = expanded_query[id];
        //my_rotate_internal(*context, temp_ct, i, *galois_keys, MemoryManager::GetPool(), NUM_EXPANSION_THREAD);
        my_rotate_internal(*context, temp_ct, i, *galois_keys, column_pools[id], NUM_EXPANSION_THREAD);
        //evaluator->add_inplace(expanded_query[id], temp_ct);
        my_add_inplace(*context, expanded_query[id], temp_ct);
    }
    return NULL;
}

void *process_rows(void *arg)
{
    int id = *((int *)arg);
    chrono::high_resolution_clock::time_point time_start, time_end, total_start, total_end;

    total_start = chrono::high_resolution_clock::now();
    Ciphertext column_results[NUM_COL_THIS_WORKER];
    vector<column_thread_arg> column_args;
    vector<mult_thread_arg> mult_args;

    for(int i = 0; i < NUM_COL_THREAD; i++) {
        column_args.push_back(column_thread_arg(i, id, column_results));
    }
    for(int i = 0; i < NUM_COL_THIS_WORKER; i++) {
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
            if (pthread_create(&(col_process_thread[i]), NULL, process_columns, (void *)&(column_args[i])))
            {
                printf("Error creating column processing thread");
            }
        }

        for (int i = 0; i < NUM_COL_THREAD; i++)
        {
            pthread_join(col_process_thread[i], NULL);
        }
        // time_end = chrono::high_resolution_clock::now();
        //col_processing_time[row_idx] = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();

        // time_start = chrono::high_resolution_clock::now();

        for(int diff = 2; diff <= NUM_COL_THIS_WORKER; diff*=2) {

            // #pragma omp parallel for
            // for(int j = 0; j < NUM_COL_THIS_WORKER; j+=diff) {
            //     my_bfv_multiply(*context, column_results[j], column_results[j + (diff/2)]);
            //     my_relinearize_internal(*context, column_results[j] ,*relin_keys, 2, MemoryManager::GetPool());
            // }

            for(int i = 0; i < mult_args.size(); i++) {
                mult_args[i].diff = diff;
            }
            for (int i = 0; i < NUM_COL_THIS_WORKER; i+=diff)
            {
                if (pthread_create(&(col_mult_thread[i]), NULL, multiply_columns, (void *)&(mult_args[i])))
                {
                    printf("Error creating column processing thread");
                }
            }

            for (int i = 0; i < NUM_COL_THIS_WORKER; i+=diff)
            {
                pthread_join(col_mult_thread[i], NULL);
            }

        }

        // time_end = chrono::high_resolution_clock::now();
        if(!GROUP_LEADER) {
            vector<uint64_t> ct(SMALL_COEFF_COUNT, 0);
            std::copy(column_results[0].data(), column_results[0].data() + SMALL_COEFF_COUNT, ct.begin());
            group_client[0]->async_call("sendGroupResponse", ct);
            return NULL;
        }

        int group_response_count = 0;
        while(group_response_count < (GROUP_SIZE - 1)) {
            sem_wait(&group_response_sem);
            pthread_mutex_lock(&group_response_lock);
            Ciphertext ct = group_response_queue.front();
            group_response_queue.pop();
            pthread_mutex_unlock(&group_response_lock);
            my_bfv_multiply(*context, column_results[0], ct, column_pools[0], TOTAL_MACHINE_THREAD);
            my_relinearize_internal(*context, column_results[0] ,*relin_keys, 2, column_pools[0], TOTAL_MACHINE_THREAD); 
            ++group_response_count;
        }

        time_start = chrono::high_resolution_clock::now();
        // evaluator->transform_to_ntt_inplace(column_results[0]);
        // evaluator->multiply_plain_inplace(column_results[0], mult_pt[row_idx]);
        Ciphertext temp_ct = column_results[0];
        // evaluator->rotate_columns(column_results[0], *galois_keys, temp_ct);
        my_conjugate_internal(*context, temp_ct, *galois_keys, column_pools[0], TOTAL_MACHINE_THREAD/NUM_ROW_THREAD);

        my_bfv_multiply(*context, column_results[0], temp_ct, column_pools[0], TOTAL_MACHINE_THREAD/ NUM_ROW_THREAD);
        my_relinearize_internal(*context, column_results[0] ,*relin_keys, 2, MemoryManager::GetPool(), TOTAL_MACHINE_THREAD / NUM_ROW_THREAD);

        my_transform_to_ntt_inplace(*context, column_results[0], (TOTAL_MACHINE_THREAD / NUM_ROW_THREAD));
        // my_multiply_plain_ntt(*context, column_results[0], mult_pt[row_idx], (TOTAL_MACHINE_THREAD / NUM_ROW_THREAD));

        row_result[row_idx] = column_results[0];
 
        time_end = chrono::high_resolution_clock::now();
        //row_agg_time[row_idx] = (chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count();
    }
    total_end = chrono::high_resolution_clock::now();

    //row_thread_processing_time_arr[id] = (chrono::duration_cast<chrono::microseconds>(total_end - total_start)).count();

}

void *multiply_columns(void *arg) {

    mult_thread_arg mult_arg = *((mult_thread_arg *)arg);
    Ciphertext *column_results = mult_arg.column_result;
    int id = mult_arg.id;
    int diff = mult_arg.diff;
    int num_threads = TOTAL_MACHINE_THREAD / (NUM_COL_THIS_WORKER / diff);

    my_bfv_multiply(*context, column_results[id], column_results[id + (diff/2)], column_pools[id], num_threads);

    my_relinearize_internal(*context, column_results[id] ,*relin_keys, 2, column_pools[id], num_threads); 
}


void *process_columns(void *arg) {

    column_thread_arg col_arg = *((column_thread_arg *)arg);
    chrono::high_resolution_clock::time_point time_start, time_end, total_start, total_end;
    vector<unsigned long> exp_time;
    unsigned long exponent_time = 0;
    int num_col_per_thread = NUM_COL_THIS_WORKER / NUM_COL_THREAD;
    int start_idx = num_col_per_thread * col_arg.col_id;
    int end_idx = start_idx + num_col_per_thread;
    time_start = chrono::high_resolution_clock::now();
    for(int i = start_idx; i < end_idx; i++) {
        Ciphertext sub;
        Ciphertext prod;
        evaluator->sub_plain(expanded_query[i], db[col_arg.row_idx][i], sub);


        for (int k = 0; k < 16; k++){
            //time_start = chrono::high_resolution_clock::now();

            //evaluator->square_inplace(sub);
            my_bfv_square(*context, sub, column_pools[i], NUM_EXPONENT_THREAD);
            //time_end = chrono::high_resolution_clock::now();
            //square_time.push_back((chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count());

            //evaluator->relinearize_inplace(sub, *relin_keys);
            my_relinearize_internal(*context, sub, *relin_keys, 2, column_pools[i], NUM_EXPONENT_THREAD);
        }
        //auto time_end = chrono::high_resolution_clock::now();
        //exp_time.push_back((chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count());

        for(int k = 0; k < MOD_SWITCH_COUNT; k++) {
            my_mod_switch_scale_to_next(*context, sub, sub, column_pools[i], NUM_EXPONENT_THREAD);
        }
        evaluator->sub(*one_ct, sub, (col_arg.column_result)[i]);

    }
    time_end = chrono::high_resolution_clock::now(); 
}

void sendGroupResponse(vector<uint64_t> ct) {
    assert(ct.size() == SMALL_COEFF_COUNT);
    Ciphertext response = *one_ct;
    std::copy(ct.begin(), ct.end(), response.data());
    pthread_mutex_lock(&group_response_lock);
    group_response_queue.push(response);
    pthread_mutex_unlock(&group_response_lock);
    sem_post(&group_response_sem);
    return;
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
        my_rotate_internal(*context, right_sum, -mid, *galois_keys, column_pools[0], TOTAL_MACHINE_THREAD);
        //evaluator->rotate_rows_inplace(right_sum, -mid, *galois_keys);
        //evaluator->add_inplace(left_sum, right_sum);
        my_add_inplace(*context, left_sum, right_sum);
        return left_sum;
    }
    else
    {

        seal::Ciphertext column_sum = query[0];
        seal::Ciphertext temp_ct;
        //evaluator->multiply_plain(query[0], pir_encoded_db[pir_num_query_ciphertext * start], column_sum);
        my_multiply_plain_ntt(*context, column_sum, pir_encoded_db[pir_num_query_ciphertext * start], TOTAL_MACHINE_THREAD);

        for (int j = 1; j < pir_num_query_ciphertext; j++)
        {
            temp_ct = query[j];
            my_multiply_plain_ntt(*context, temp_ct, pir_encoded_db[pir_num_query_ciphertext * start + j], TOTAL_MACHINE_THREAD);

            // evaluator->multiply_plain(query[j], pir_encoded_db[pir_num_query_ciphertext * start + j], temp_ct);
            // evaluator->add_inplace(column_sum, temp_ct);
            my_add_inplace(*context, column_sum, temp_ct);
        }
        my_transform_from_ntt_inplace(*context, column_sum, TOTAL_MACHINE_THREAD);
        //evaluator->transform_from_ntt_inplace(column_sum);
        return column_sum;
    }
}

void init_pir_params()
{
    pir_plain_bit_count = PLAIN_BIT;
    pir_num_obj = ((N/2) * NUM_ROW);
    pir_num_query_ciphertext = ceil(pir_num_obj / (double)(N/2));
    pir_num_columns_per_obj = 2 * (ceil(((pir_obj_size/2) * 8) / (float)(pir_plain_bit_count)));
    pir_db_rows = ceil(pir_num_obj / (double)N) * pir_num_columns_per_obj;

    return;
}

void populate_pir_db() {

    vector<vector<uint64_t>> pir_db; //64 bit placeholder for 16 bit plaintext coefficients
    for(int i = 0; i < pir_num_obj; i++) {
        vector<uint64_t> v;
        for(int j = 0; j < (pir_obj_size / 2); j++) {  // 2 bytes each plaintxt slot
            v.push_back(((i+1)*(j+1) + i*i + (j+1)*(j+1)+10) % PLAIN_MODULUS);
        }
        pir_db.push_back(v);
    }

    set_pir_db(pir_db);

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


void preallocate_memory() {

    #pragma omp parallel for num_threads(NUM_COL_THIS_WORKER)
    for(int i = 0; i < NUM_COL_THIS_WORKER; i++) {
        SEAL_ALLOCATE_GET_POLY_ITER(temp, 3, N, 13, column_pools[i]);
        SEAL_ALLOCATE_GET_POLY_ITER(temp2, 3, N, 14, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp4, (8192*4), 15, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp5, (8192*4), 14, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp6, (8192*4), 13, column_pools[i]);
        SEAL_ALLOCATE_GET_COEFF_ITER(r_m_tilde,N, column_pools[i]);
        SEAL_ALLOCATE_GET_COEFF_ITER(temp7, N, column_pools[i]);

        SEAL_ALLOCATE_GET_STRIDE_ITER(temp8, uint64_t, N, 13, column_pools[i]);
        SEAL_ALLOCATE_GET_STRIDE_ITER(temp9, uint64_t, N, 5, column_pools[i]);

        SEAL_ALLOCATE_GET_RNS_ITER(temp10, N, 27, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp11, N, 26, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp12, N, 24, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp13, N, 22, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp14, N, 20, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp15, N, 18, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp16, N, 16, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp17, N, 12, column_pools[i]);


        auto t_poly_prod(allocate_zero_poly_array(2, N, 14, column_pools[i]));
        auto t_poly_lazy(allocate_zero_poly_array(2, N, 2, column_pools[i]));
        auto t_poly_lazy2(allocate_zero_poly_array(2, N, 2, column_pools[i]));
        auto t_poly_lazy3(allocate_zero_poly_array(2, N, 2, column_pools[i]));

        SEAL_ALLOCATE_GET_RNS_ITER(temp18, N, 7, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp19, N, 7, column_pools[i]);

        SEAL_ALLOCATE_GET_POLY_ITER(encrypted1_q, 2, N, 5, column_pools[i]);
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted1_Bsk, 2, N, 6, column_pools[i]);
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted1_q2, 2, N, 5, column_pools[i]);

        SEAL_ALLOCATE_GET_RNS_ITER(temp_q_Bsk, N, 11, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp_q_Bsk2, N, 6, column_pools[i]);
        SEAL_ALLOCATE_GET_RNS_ITER(temp20, N, 1, column_pools[i]);

        SEAL_ALLOCATE_GET_POLY_ITER(temp21, 1, N, 1, column_pools[i]);

        auto t_poly_prod2(allocate_zero_poly_array(2, N, 2, column_pools[i]));

    }
}
