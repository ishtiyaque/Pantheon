#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>
#include <sstream> 
#include<cmath>
#include<ctime>
#include<stack>
#include <unistd.h> 
#include<future>
#include<cassert>
#include<set>
#include<iostream>
#include "seal/seal.h"
//#include "md5.h"

#include "utils.hpp"

#include<rpc/client.h>
#include<rpc/server.h>

using namespace std::chrono;
using namespace std;
using namespace seal;
    
#define N 32768
#define MOD_SWITCH_COUNT 9

#define MASTER_PORT 4000
#define CLIENT_PORT 2000
#define WORKER_PORT 3000

#define PLAIN_BIT 16
#define PLAIN_MODULUS 65537
vector<int> CT_PRIMES({60,60,60,60,60,60,60,60,60,60,60,60,60});

#define LARGE_COEFF_COUNT (((CT_PRIMES.size() - 1) * (N) * 2))
#define SMALL_COEFF_COUNT (((CT_PRIMES.size() - 1 - MOD_SWITCH_COUNT) * (N) * 2))

typedef  std::vector<seal::Ciphertext> PIRQuery;
typedef  seal::Ciphertext PIRResponse;
