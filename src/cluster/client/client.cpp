#include<globals.h>

using namespace std::chrono;
using namespace std;
using namespace seal;

string CLIENT_IP = "";
string MASTER_IP = "";

int NUM_COL = 0;
uint64_t  latency;

uint64_t query_send_timestamp;
uint64_t response_rcv_timestamp;

uint64_t query_gen_time;
uint64_t response_decode_time;
uint64_t round_trip_time;

uint32_t pir_obj_size = 0;
int pir_num_columns_per_obj;

size_t query_size, response_size;

clock_t total_cpu_start_time, total_cpu_stop_time;

void sendResponse(int id,string response);
void printReport();


int main(int argc, char *argv[]) {

    uint64_t number_of_items = 0;

    int option;
    const char *optstring = "n:p:c:q:s:";
    while ((option = getopt(argc, argv, optstring)) != -1) {
	switch (option) {

        case 'n':
            number_of_items = stoi(optarg);
            break;
        case 'p':
            MASTER_IP = string(optarg);
            break;
        case 'q':
            CLIENT_IP = string(optarg);
            break;

        case 'c':
            NUM_COL = stoi(optarg);
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
    if (MASTER_IP.size() < 7){cout << "Missing -p\n"; return 0;}
    if(CLIENT_IP.size() < 7) {cout<<"Missing -q\n";return 0;}
    if(!NUM_COL) {cout<<"Missing -c\n";return 0;} 
    if(!pir_obj_size) {cout<<"Missing -s\n";return 0;} 

    pir_num_columns_per_obj = 2 * (ceil(((pir_obj_size/2) * 8) / (float)(PLAIN_BIT)));

    rpc::client query_client(MASTER_IP, MASTER_PORT);

    chrono::high_resolution_clock::time_point time_start, time_end, total_start, total_end;

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(N);
    parms.set_coeff_modulus(CoeffModulus::Create(N, CT_PRIMES));
    parms.set_plain_modulus(PLAIN_MODULUS);

    SEALContext context(parms);
    auto pid = context.first_parms_id();

    uint64_t plain_modulus = parms.plain_modulus().value();
    // cout << parms.plain_modulus().value() << endl;
    // cout << parms.coeff_modulus().at(0).value() << endl;
    // cout << parms.coeff_modulus().at(1).value() << endl;
    // print_parameters(context);
    //cout << endl;
    // cout << "Maximum size of co-efficient modulus: " << CoeffModulus::MaxBitCount(N) << endl;
    //auto qualifiers = context->first_context_data()->qualifiers();
    //cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;
    //cout << "Fast plain lift enabled: " << qualifiers.using_fast_plain_lift << endl;

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();  

    set<int> rotation_steps;  
    rotation_steps.insert(0);

    for (int i = N / (2 * NUM_COL); i < N / 2; i *= 2) {
        rotation_steps.insert(i);
    }

    for (int i = 1; i <= (pir_num_columns_per_obj / 2); i *= 2)
    {
        rotation_steps.insert(-i);
    }
    // gal_keys = keygen->galois_keys_local(steps);

    Serializable <GaloisKeys> galois_keys = keygen.create_galois_keys(vector<int>(rotation_steps.begin(), rotation_steps.end()));
    //Serializable <GaloisKeys> galois_keys = keygen.create_galois_keys();

    Serializable <RelinKeys> relin_keys = keygen.create_relin_keys();

    Encryptor encryptor(context, secret_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();
    size_t row_size = slot_count / 2;

    vector<uint64_t> one_mat;
    Plaintext one_pt;
    for (int i = 0; i < N; i++)
    {
        one_mat.push_back(1);
    }
    batch_encoder.encode(one_mat, one_pt);
    Serializable<Ciphertext> one_ct = encryptor.encrypt_symmetric(one_pt);
    stringstream oss;
    one_ct.save(oss);

    std::stringstream gkss, rkss, qss;
    string serialized_gal_key, serialized_relin_key, serialized_one_ct ;
    galois_keys.save(gkss);
    relin_keys.save(rkss);
    serialized_gal_key = gkss.str();
    serialized_relin_key = rkss.str();
    serialized_one_ct = oss.str();

    vector<uint64_t> query_mat(N, 0), response_mat(N, 0);
    Plaintext query_pt, response_pt;
    Ciphertext response;  
    response.load(context, oss);
    for (int k = 0; k < MOD_SWITCH_COUNT; k++)
    {
        evaluator.mod_switch_to_next_inplace(response);
    }

    query_client.call("sendKeys", serialized_gal_key, serialized_relin_key, serialized_one_ct);
    sleep(50);
 
    total_cpu_start_time = clock();

    total_start = chrono::high_resolution_clock::now();

    {
        int select_row = 0;
        int select_idx = 5;
        for(int i = 0; i < N/(NUM_COL*2);i++) {
            for(int j = 0; j < NUM_COL; j++) {
                query_mat[j*(N/(NUM_COL*2)) + i] = (j + select_row * NUM_COL + select_idx) % PLAIN_MODULUS;
                query_mat[j*(N/(NUM_COL*2)) + i + N/2] = (j + select_row * NUM_COL + select_idx) % PLAIN_MODULUS;
            }
        }
    }

    batch_encoder.encode(query_mat, query_pt);

    Serializable <Ciphertext> ser_query = encryptor.encrypt_symmetric(query_pt);
    ser_query.save(qss);
    string serialized_query = qss.str();
    time_end = chrono::high_resolution_clock::now();
    query_gen_time = chrono::duration_cast<chrono::microseconds>(time_end - total_start).count();
    query_send_timestamp = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    vector<uint64_t> ct = query_client.call("sendQuery", serialized_query).as<vector<uint64_t>>();
    response_rcv_timestamp = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();

    assert(ct.size() == SMALL_COEFF_COUNT );
    std::copy(ct.begin(), ct.end(), response.data());
    decryptor.decrypt(response, response_pt);
    batch_encoder.decode(response_pt, response_mat);
 
    total_end = chrono::high_resolution_clock::now();
    total_cpu_stop_time = clock();

    response_decode_time = chrono::duration_cast<chrono::microseconds>(total_end.time_since_epoch()).count() - response_rcv_timestamp;
    latency = chrono::duration_cast<chrono::microseconds>(total_end - total_start).count();
    query_size = serialized_query.size();
    response_size = (ct.size() * sizeof(uint64_t));
    cout<<"Response noise budget"<<endl<<decryptor.invariant_noise_budget(response)<<endl<<endl;
    printReport();

    return 0;
}



void printReport()
{
    cout<<"latency (us)"<<endl<<latency<<endl<<endl;

    cout << "query send timestamp:" << endl;  
    cout<<query_send_timestamp<<endl;

    cout << "\nresponse_received timestamp: " << endl;
    cout << response_rcv_timestamp<<endl;

    cout<<"\nTotal CPU time (sec): "<<endl<<(((float)total_cpu_stop_time-total_cpu_start_time)/CLOCKS_PER_SEC)<<endl<<endl;
    
    cout<<"query gen time "<<endl;
    cout<<query_gen_time<<endl;

    cout<<"\nround trip time "<<endl;
    cout<<response_rcv_timestamp - query_send_timestamp<<endl;

    cout<<"\nresponse decode time "<<endl;
    cout<<response_decode_time<<endl;

    cout<<"\nquery size (byte)"<<endl<<query_size<<endl;
    cout<<"\nresponse size (byte)"<<endl<<response_size<<endl;

    return;
}