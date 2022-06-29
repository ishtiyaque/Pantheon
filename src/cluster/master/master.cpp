#include<globals.h>

pthread_mutex_t request_received_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t request_received_cond = PTHREAD_COND_INITIALIZER;

pthread_mutex_t gal_key_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t gal_key_cond = PTHREAD_COND_INITIALIZER;

pthread_mutex_t worker_response_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t received_all_response = PTHREAD_COND_INITIALIZER;

pthread_mutex_t sent_client_response_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t sent_client_response_cond = PTHREAD_COND_INITIALIZER;

SEALContext *context;
Evaluator *evaluator;

vector<uint64_t> sendQuery(string _serialized_query);
void sendKeys( string _serialized_gal_key, string _serialized_relin_key, string _serialized_one_ct);
void *sendToWorker(void *arg);
void prepareSendThread();
void readWorkerData();
void printReport();

string serialized_gal_key, serialized_relin_key;
string serialized_one_ct;
string serialized_query;
string MASTER_IP;
int NUM_TOTAL_WORKER = 0;
int GROUP_SIZE = 0;
int NUM_GROUP;
int worker_response_count = 0;

uint64_t query_rcv_timestamp, response_send_timestamp;

vector<Ciphertext> worker_response;

clock_t total_cpu_start_time, total_cpu_stop_time;

rpc::client **worker_client;

pthread_t *send_threads;
int *thread_id;
string *worker_ip;
int main(int argc, char **argv) {

    int option;
    const char *optstring = "p:w:g:";
    while ((option = getopt(argc, argv, optstring)) != -1) {
	switch (option) {
       case 'p':
            MASTER_IP = string(optarg);
            break;
        case 'w':
            NUM_TOTAL_WORKER = stoi(optarg);
            break;
        case 'g':
            GROUP_SIZE = stoi(optarg);
            break;  

        case '?':
            cout<<"error optopt: "<<optopt<<endl;
            cout<<"error opterr: "<<opterr<<endl;
            return 1;
	    }
    }
    if(MASTER_IP.size() < 7) {cout<<"Missing -p\n";return 0;}
    if(!NUM_TOTAL_WORKER) {cout<<"Missing -w\n";return 0;}
    if(!GROUP_SIZE) {cout<<"Missing -g\n"; return 0;}

    NUM_GROUP = (NUM_TOTAL_WORKER / GROUP_SIZE);

    rpc::server query_server(string(MASTER_IP), MASTER_PORT);
    query_server.bind("sendKeys", sendKeys);
    query_server.bind("sendQuery", sendQuery);
    readWorkerData();

    srand(time(NULL));

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(N);

    //parms.set_coeff_modulus(CoeffModulus::BFVDefault(N));

    parms.set_coeff_modulus(CoeffModulus::Create(N, CT_PRIMES));
    

    //parms.set_plain_modulus(PlainModulus::Batching(N, 15+1));
    parms.set_plain_modulus(PLAIN_MODULUS);

    //parms.set_plain_modulus(17);
    context = new SEALContext(parms);
    auto pid = context->first_parms_id();
    evaluator = new Evaluator(*context);

    query_server.async_run(1);

    pthread_mutex_lock(&gal_key_lock);
    pthread_cond_wait(&gal_key_cond, &gal_key_lock);
    pthread_mutex_unlock(&gal_key_lock);

    Ciphertext one_ct;
    stringstream oss(serialized_one_ct);
    one_ct.load(*context, oss);
    for (int k = 0; k < MOD_SWITCH_COUNT; k++)
    {
        evaluator->mod_switch_to_next_inplace(one_ct);
    }
    for(int i = 0; i < NUM_GROUP; i++) {
        worker_response.push_back(one_ct);
    }

    prepareSendThread();   
    total_cpu_start_time = clock();

    pthread_mutex_lock(&sent_client_response_lock);
    pthread_cond_wait(&sent_client_response_cond, &sent_client_response_lock);
    pthread_mutex_unlock(&sent_client_response_lock);

    sleep(3);
    total_cpu_stop_time = clock();

    printReport();
    return 0;
}


vector<uint64_t> sendQuery(string _serialized_query) {
    // cout<<"received query"<<endl;
    query_rcv_timestamp = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    serialized_query = _serialized_query;
    Ciphertext result;
    vector<uint64_t> response(SMALL_COEFF_COUNT, 0);
    
    pthread_mutex_lock(&request_received_lock);
    pthread_cond_broadcast(&request_received_cond);
    pthread_mutex_unlock(&request_received_lock);

    pthread_mutex_lock(&worker_response_lock);
    pthread_cond_wait(&received_all_response, &worker_response_lock);
    pthread_mutex_unlock(&worker_response_lock);
    result = worker_response[0];
    for(int i = 1; i < worker_response.size(); i++) {
        // evaluator->add_inplace(result, worker_response[i]);  // TO DO : parallelize the addition
        my_add_inplace(*context, result, worker_response[i]);
    }

    std::copy(result.data(), result.data() + SMALL_COEFF_COUNT, response.begin()); 

    pthread_mutex_lock(&sent_client_response_lock);
    pthread_cond_broadcast(&sent_client_response_cond);
    pthread_mutex_unlock(&sent_client_response_lock);

    response_send_timestamp = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    return response;
}

void sendKeys( string _serialized_gal_key, string _serialized_relin_key, string _serialized_one_ct) {
    // cout<<"received keys"<<endl;
    serialized_gal_key = _serialized_gal_key;
    serialized_relin_key = _serialized_relin_key;
    serialized_one_ct = _serialized_one_ct;
    pthread_mutex_lock(&gal_key_lock);
    pthread_cond_broadcast(&gal_key_cond);
    pthread_mutex_unlock(&gal_key_lock);
    return;
}

void prepareSendThread() {

    thread_id = new int[NUM_TOTAL_WORKER];
    send_threads = new pthread_t[NUM_TOTAL_WORKER];

    for (int i = 0; i < NUM_TOTAL_WORKER; i++)
    {
        thread_id[i] = i;
        if (pthread_create(&send_threads[i], NULL, sendToWorker, (void *)(&thread_id[i])))
        {            
            printf("Error creating thread\n");
            exit(1);
        }
        if(i%4 == 0) sleep(1);
    }
    // for (int i = 0; i < NUM_TOTAL_WORKER; i++)
    // {
    //     pthread_join(send_threads[i], NULL);
    // }

    return;
}

void *sendToWorker(void *arg) {
    int my_id = * ((int *)arg);
    rpc::client worker_client(worker_ip[my_id], WORKER_PORT + my_id);
    //cout<<"connected to worker at "<<worker_ip[my_id]<<":"<<MASTER_WORKER_PORT + my_id<<endl;
    worker_client.call("sendKeys", serialized_gal_key, serialized_relin_key, serialized_one_ct);
    if((my_id % GROUP_SIZE) != 0) {
        return NULL;
    }
    pthread_mutex_lock(&request_received_lock);
    pthread_cond_wait(&request_received_cond, &request_received_lock);
    pthread_mutex_unlock(&request_received_lock);
    //query_send_start_timestamp[my_id] = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    vector<uint64_t> ct = worker_client.call("sendQuery", serialized_query).as<vector<uint64_t>>();
    std::copy(ct.begin(), ct.end(), worker_response[my_id / GROUP_SIZE].data());
    pthread_mutex_lock(&worker_response_lock);
    //worker_response.push_back(response);
    worker_response_count++;
    if(worker_response_count == NUM_GROUP) {
        pthread_cond_broadcast(&received_all_response);
    }
    pthread_mutex_unlock(&worker_response_lock);
    //query_send_end_timestamp[my_id] = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    return NULL;
}

void readWorkerData() {
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

    cout<<"\nTime at master"<<endl<<(response_send_timestamp - query_rcv_timestamp)<<endl;

    cout<<"\nTotal CPU time:(sec)"<<endl<<((float)total_cpu_stop_time-total_cpu_start_time)/CLOCKS_PER_SEC<<endl;

    return ;
}

