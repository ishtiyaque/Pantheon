# Pantheon: Private Retrieval from Public Key-Value Store

Consider a cloud server that owns a key-value store and provides a private query service to its clients.
The key-value store in this case is _public_ and a client cannot encrypt or modify it to ensure privacy.
Pantheon is a system that allows its clients to retrieve the value corresponding to a key from such a key-value store privately, i.e. 
without letting anyone else in the world (even the server) to know which key the client is interested in.

## Environment setup

The codebase is tested on Ubuntu 18.04 and 20.04.
To install all the dependencies, first clone this repositry:

    git clone https://github.com/ishtiyaque/Pantheon.git

Then run the following commands:

    cd Pantheon
    ./setup.sh
    
## Single-server deployment

The code for deploying Pantheon server on a single machine is in `Pantheon/src/single_machine/`

### Build
To build the code, run the following:

    cd Pantheon/src/single_machine/
    cmake .
    make
    
    
### Run

To run Pantheon's single-server implementation, run the following commands:

    cd Pantheon/src/single_machine/bin/
    ./Pantheon -n <number of items> -k <key_size> -s <value_size>
    
The options are explained below:

    -n The number of key-value tuples in the key-value store.
    -k The size of each key in bits.
    -s The size of each value in Bytes.
    
For example, to perform a private query over a key-value store containing `32768` tuples where each key is `64-bits` and each value is `256 Bytes`, run the follwoing command:

    ./Pantheon -n 32768 -k 64 -s 256

