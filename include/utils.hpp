#pragma once

#include <vector>
#include <memory>
#include <map>
#include "seal/context.h"
#include "seal/relinkeys.h"
//#include "seal/smallmodulus.h"
#include "seal/memorymanager.h"
#include "seal/ciphertext.h"
#include "seal/plaintext.h"
#include "seal/galoiskeys.h"
#include "seal/util/pointer.h"
#include "seal/secretkey.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/common.h"
#include "seal/evaluator.h"
#include "seal/util/common.h"
#include "seal/util/galois.h"
#include "seal/util/numth.h"
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include "seal/util/defines.h"

#include "seal/util/uintarith.h"
#include <algorithm>
#include <cmath>
#include "seal/util/scalingvariant.h"
#include <functional>
#include<iostream>
#include <chrono>

#include "omp.h"

using namespace seal;
using namespace seal::util;

using namespace std;

//int NUM_OMP_THREAD = 4;
void my_add_inplace(SEALContext &context_, Ciphertext &encrypted1, Ciphertext &encrypted2);
void my_bfv_square(SEALContext &context_, Ciphertext &encrypted, MemoryPoolHandle pool, int num_threads);
void my_fastbconv_m_tilde(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads);
void my_fast_convert_array(const Pointer<BaseConverter> &conv, ConstRNSIter in, RNSIter out, MemoryPoolHandle pool, int num_threads);
inline void my_multiply_poly_scalar_coeffmod(ConstRNSIter poly, std::size_t coeff_modulus_size, std::uint64_t scalar, 
                                                ConstModulusIter modulus,RNSIter result, int num_threads);
void my_sm_mrq(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads);
void my_fast_floor(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads);
void my_fastbconv_sk(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads);

void my_relinearize_internal(SEALContext &context_, Ciphertext &encrypted, const RelinKeys &relin_keys, size_t destination_size, MemoryPoolHandle pool, int num_threads);
// void my_switch_key_inplace( SEALContext &context_,
//         Ciphertext &encrypted, ConstRNSIter target_iter, const KSwitchKeys &kswitch_keys, size_t kswitch_keys_index,
//         MemoryPoolHandle pool);
void my_bfv_multiply(SEALContext &context_, Ciphertext &encrypted1, const Ciphertext &encrypted2, MemoryPoolHandle pool, int num_threads);
void my_mod_switch_scale_to_next( SEALContext &context_, Ciphertext &encrypted, Ciphertext &destination, MemoryPoolHandle pool, int num_threads);
void my_transform_to_ntt_inplace(SEALContext &context_, Ciphertext &encrypted, int num_threads);
void my_transform_from_ntt_inplace(SEALContext &context_, Ciphertext &encrypted_ntt, int num_threads);
void my_multiply_plain_ntt(SEALContext &context_, Ciphertext &encrypted_ntt, const Plaintext &plain_ntt, int num_threads);
void my_rotate_internal(SEALContext context_, Ciphertext &encrypted, int steps, const GaloisKeys &galois_keys, MemoryPoolHandle pool, int num_threads);
void my_conjugate_internal(SEALContext context_, Ciphertext &encrypted, const GaloisKeys &galois_keys, MemoryPoolHandle pool, int num_threads) ;
void my_apply_galois_inplace(SEALContext context_, Ciphertext &encrypted, uint32_t galois_elt, const GaloisKeys &galois_keys, MemoryPoolHandle pool, int num_threads);
void my_switch_key_inplace(
    SEALContext &context_, Ciphertext &encrypted, ConstRNSIter target_iter, const KSwitchKeys &kswitch_keys, size_t kswitch_keys_index,
    MemoryPoolHandle pool, int num_threads);

inline void my_inverse_ntt_negacyclic_harvey(PolyIter operand, std::size_t size, ConstNTTTablesIter tables)
{
    //SEAL_ITERATE(operand, size, [&](auto I) { inverse_ntt_negacyclic_harvey(I, operand.coeff_modulus_size(), tables);});
    for(int i = 0; i < size;i++) {
        for(int j = 0; j < operand.coeff_modulus_size();j++) {
            inverse_ntt_negacyclic_harvey_lazy(operand[i][j], tables[j]);
        }
    }
}


void my_add_inplace(SEALContext &context_, Ciphertext &encrypted1, Ciphertext &encrypted2)
    {
        // Verify parameters.
        if (encrypted1.parms_id() != encrypted2.parms_id())
        {
            throw invalid_argument("encrypted1 and encrypted2 parameter mismatch");
        }
        if (encrypted1.is_ntt_form() != encrypted2.is_ntt_form())
        {
            throw invalid_argument("NTT form mismatch");
        }
        if(encrypted1.size() != encrypted2.size()) {
            throw invalid_argument("encrypted1 and encrypted2 size mismatch");
        }

        // Extract encryption parameters.
        auto &context_data = *context_.get_context_data(encrypted1.parms_id());
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t encrypted_size = encrypted1.size();
        auto iter1 = PolyIter(encrypted1);
        auto iter2 = PolyIter(encrypted2);

        #pragma omp parallel for collapse(2)
        for(int i = 0; i < encrypted_size; i++) {
            for(int j = 0; j < coeff_modulus_size; j++) {
                add_poly_coeffmod(iter1[i][j], iter2[i][j], coeff_count, coeff_modulus[j], iter1[i][j]);
            }
        }
    }

inline void my_add_poly_coeffmod(
    ConstRNSIter operand1, ConstRNSIter operand2, std::size_t coeff_modulus_size, ConstModulusIter modulus,
    RNSIter result, int num_threads)
{
    auto poly_modulus_degree = result.poly_modulus_degree();

    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0;i< coeff_modulus_size;i++) {
        add_poly_coeffmod(operand1[i], operand2[i], poly_modulus_degree, modulus[i], result[i]);
    }
}

inline void my_dyadic_product_coeffmod(
    ConstRNSIter operand1, ConstRNSIter operand2, std::size_t coeff_modulus_size, ConstModulusIter modulus,
    RNSIter result, int num_threads)
{
    auto poly_modulus_degree = result.poly_modulus_degree();
    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < coeff_modulus_size; i++) {
        dyadic_product_coeffmod(operand1[i], operand2[i], poly_modulus_degree, modulus[i], result[i]);
    }
}

inline void my_ntt_negacyclic_harvey_lazy(
    RNSIter operand, std::size_t coeff_modulus_size, ConstNTTTablesIter tables, int num_threads)
{
    #pragma omp paralle for num_threads(num_threads)
    for(int i = 0; i < coeff_modulus_size ;i++) {
        ntt_negacyclic_harvey_lazy(operand[i], tables[i]);
    }
}

void my_bfv_square(SEALContext &context_, Ciphertext &encrypted, MemoryPoolHandle pool, int num_threads) 

    {
        omp_set_num_threads(num_threads);
        chrono::high_resolution_clock::time_point time_start, time_end, total_start, total_end;
time_start = chrono::high_resolution_clock::now();

        if (encrypted.is_ntt_form())
        {
            throw invalid_argument("encrypted cannot be in NTT form");
        }

        // Extract encryption parameters.
        auto &context_data = *context_.get_context_data(encrypted.parms_id());
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree(); // N
        size_t base_q_size = parms.coeff_modulus().size(); // 13
        size_t encrypted_size = encrypted.size(); // 2
        uint64_t plain_modulus = parms.plain_modulus().value();

        double new_scale = encrypted.scale();
        // if (!is_scale_within_bounds(new_scale, context_data))
        // {
        //     //throw invalid_argument("scale out of bounds");
        // }

        auto rns_tool = context_data.rns_tool();
        size_t base_Bsk_size = rns_tool->base_Bsk()->size(); // 14
        size_t base_Bsk_m_tilde_size = rns_tool->base_Bsk_m_tilde()->size(); // 15

        // Optimization implemented currently only for size 2 ciphertexts
        // if (encrypted_size != 2)
        // {
        //     bfv_multiply(encrypted, encrypted, move(pool));
        //     return;
        // }

        // Determine destination.size()
        //size_t dest_size = sub_safe(add_safe(encrypted_size, encrypted_size), size_t(1));
        size_t dest_size = encrypted_size + encrypted_size - 1;


        // Size check
        //if (!product_fits_in(dest_size, coeff_count, base_Bsk_m_tilde_size))
        //{
            //throw logic_error("invalid parameters");
        //}

        // Set up iterators for bases
        auto base_q = iter(parms.coeff_modulus());
        auto base_Bsk = iter(rns_tool->base_Bsk()->base());


        // Set up iterators for NTT tables
        auto base_q_ntt_tables = iter(context_data.small_ntt_tables());
        auto base_Bsk_ntt_tables = iter(rns_tool->base_Bsk_ntt_tables());

        // Microsoft SEAL uses BEHZ-style RNS multiplication. For details, see Evaluator::bfv_multiply. This function
        // uses additionally Karatsuba multiplication to reduce the complexity of squaring a size-2 ciphertext, but the
        // steps are otherwise the same as in Evaluator::bfv_multiply.

        // Resize encrypted to destination size

        encrypted.resize(context_, context_data.parms_id(), dest_size);

        // Allocate space for a base q output of behz_extend_base_convert_to_ntt
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted_q, encrypted_size, coeff_count, base_q_size, pool);

        // Allocate space for a base Bsk output of behz_extend_base_convert_to_ntt
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted_Bsk, encrypted_size, coeff_count, base_Bsk_size, pool);

        // Perform BEHZ steps (1)-(3)

        auto encrypted_iter = PolyIter(encrypted);
        SEAL_ALLOCATE_GET_RNS_ITER(temp, coeff_count, base_Bsk_m_tilde_size, pool);

        SEAL_ALLOCATE_GET_RNS_ITER(temp2, (8192*4), base_q_size, pool); // Need to replace the hard coded value

        SEAL_ALLOCATE_ZERO_GET_POLY_ITER(temp_dest_q, dest_size, coeff_count, base_q_size, pool);
        SEAL_ALLOCATE_ZERO_GET_POLY_ITER(temp_dest_Bsk, dest_size, coeff_count, base_Bsk_size, pool);

time_end = chrono::high_resolution_clock::now();
//cout<<"Pre time: "<<(chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count()<<endl;

        time_start = chrono::high_resolution_clock::now();

        //SEAL_ITERATE(iter(encrypted, encrypted_q, encrypted_Bsk), encrypted_size, behz_extend_base_convert_to_ntt);

        //#pragma omp parallel for    // Incorrect result
        for(int j = 0; j < encrypted_size; j++) {

            #pragma omp parallel for
            for(int i = 0; i < base_q_size;i++) { //3.7 4.6 ms
                set_uint(encrypted_iter[j][i], coeff_count, encrypted_q[j][i]);
                ntt_negacyclic_harvey_lazy(encrypted_q[j][i], base_q_ntt_tables[i]);
                multiply_poly_scalar_coeffmod(encrypted_iter[j][i], temp2.poly_modulus_degree(), rns_tool->m_tilde().value(), rns_tool->base_q()->base()[i], temp2[i]);
            }

            my_fast_convert_array(rns_tool->base_q_to_Bsk_conv(), temp2, temp, pool, num_threads);
            my_fast_convert_array(rns_tool->base_q_to_m_tilde_conv(), temp2, temp + base_Bsk_size, pool, num_threads);

            //my_sm_mrq(rns_tool, temp, encrypted_Bsk[j], pool); // 5ms

            ConstCoeffIter input_m_tilde = temp[base_Bsk_size];
            const uint64_t m_tilde_div_2 = rns_tool->m_tilde().value() >> 1;

            SEAL_ALLOCATE_GET_COEFF_ITER(r_m_tilde, rns_tool->coeff_count(), pool);
            multiply_poly_scalar_coeffmod( input_m_tilde, rns_tool->coeff_count(), rns_tool->neg_inv_prod_q_mod_m_tilde(), rns_tool->m_tilde(), r_m_tilde);

            #pragma omp parallel for    
            for(int i = 0; i < base_Bsk_size; i++) {
                MultiplyUIntModOperand prod_q_mod_Bsk_elt;
                prod_q_mod_Bsk_elt.set(rns_tool->prod_q_mod_Bsk()[i], rns_tool->base_Bsk()->base()[i]);
                SEAL_ITERATE(iter(temp[i], r_m_tilde, encrypted_Bsk[j][i]), rns_tool->coeff_count(), [&](auto J) {
                    uint64_t temp = get<1>(J);
                    if (temp >= m_tilde_div_2)
                    {
                        temp += rns_tool->base_Bsk()->base()[i].value() - rns_tool->m_tilde().value();
                    }
                    get<2>(J) = multiply_uint_mod(
                        multiply_add_uint_mod(temp, prod_q_mod_Bsk_elt, get<0>(J),rns_tool->base_Bsk()->base()[i]), rns_tool->inv_m_tilde_mod_Bsk()[i],
                        rns_tool->base_Bsk()->base()[i]);
                });
                ntt_negacyclic_harvey_lazy(encrypted_Bsk[j][i], base_Bsk_ntt_tables[i]);

            }

        }
        time_end = chrono::high_resolution_clock::now();
        //cout<<"Step 1 to 3 time: "<<(chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count()<<endl;
 
       // Allocate temporary space for the output of step (4)
        // We allocate space separately for the base q and the base Bsk components

        // Perform BEHZ step (4): dyadic Karatsuba-squaring on size-2 ciphertexts

        // This lambda function computes the size-2 ciphertext square for BFV multiplication. Since we use the BEHZ
        // approach, the multiplication of individual polynomials is done using a dyadic product where the inputs
        // are already in NTT form. The arguments of the lambda function are expected to be as follows:
        //
        // 1. a ConstPolyIter pointing to the beginning of the input ciphertext (in NTT form)
        // 3. a ConstModulusIter pointing to an array of Modulus elements for the base
        // 4. the size of the base
        // 5. a PolyIter pointing to the beginning of the output ciphertext
        auto behz_ciphertext_square = [&](ConstPolyIter in_iter, ConstModulusIter base_iter, size_t base_size,
                                          PolyIter out_iter) {
            // Compute c0^2
            //my_dyadic_product_coeffmod(in_iter[0], in_iter[0], base_size, base_iter, out_iter[0]);

            #pragma omp parallel for num_threads(num_threads)
            for(int i = 0; i < base_size; i++) {
                dyadic_product_coeffmod(in_iter[0][i], in_iter[0][i], out_iter[0].poly_modulus_degree(), base_iter[i], out_iter[0][i]);
                dyadic_product_coeffmod(in_iter[0][i], in_iter[1][i], out_iter[1].poly_modulus_degree(), base_iter[i], out_iter[1][i]);
                add_poly_coeffmod(out_iter[1][i], out_iter[1][i], out_iter[1].poly_modulus_degree(), base_iter[i], out_iter[1][i]);
                dyadic_product_coeffmod(in_iter[1][i], in_iter[1][i], out_iter[2].poly_modulus_degree(), base_iter[i], out_iter[2][i]);
            }
        };

        // Perform the BEHZ ciphertext square both for base q and base Bsk
    time_start = chrono::high_resolution_clock::now();
        behz_ciphertext_square(encrypted_q, base_q, base_q_size, temp_dest_q);
        behz_ciphertext_square(encrypted_Bsk, base_Bsk, base_Bsk_size, temp_dest_Bsk);
    time_end = chrono::high_resolution_clock::now();
    //cout<<"step 4 time: "<<(chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count()<<endl;

    SEAL_ALLOCATE_GET_RNS_ITER(temp_q_Bsk, coeff_count, base_q_size + base_Bsk_size, pool);
    SEAL_ALLOCATE_GET_RNS_ITER(temp_Bsk, coeff_count, base_Bsk_size, pool);

        // Perform BEHZ step (5): transform data from NTT form
    time_start = chrono::high_resolution_clock::now();
        //#pragma omp parallel for 
        for(int i = 0; i < dest_size;i++) {
            #pragma omp parallel for 
            for(int j = 0; j < temp_dest_Bsk.coeff_modulus_size();j++) { //14
                inverse_ntt_negacyclic_harvey_lazy(temp_dest_Bsk[i][j], base_Bsk_ntt_tables[j]);
                multiply_poly_scalar_coeffmod(temp_dest_Bsk[i][j], (temp_q_Bsk  + base_q_size).poly_modulus_degree(), plain_modulus, base_Bsk[j], (temp_q_Bsk + base_q_size)[j]);
                if(j < temp_dest_q.coeff_modulus_size()) {
                    inverse_ntt_negacyclic_harvey_lazy(temp_dest_q[i][j], base_q_ntt_tables[j]);
                    multiply_poly_scalar_coeffmod(temp_dest_q[i][j], temp_q_Bsk.poly_modulus_degree(), plain_modulus, base_q[j], temp_q_Bsk[j]);
                } 

            }
            my_fast_floor(rns_tool, temp_q_Bsk, temp_Bsk, pool, num_threads);
            my_fastbconv_sk(rns_tool, temp_Bsk, encrypted_iter[i], pool, num_threads);

        }

    time_end = chrono::high_resolution_clock::now();
    //cout<<"step 5 time: "<<(chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count()<<endl;

        // Perform BEHZ steps (6)-(8)
    time_start = chrono::high_resolution_clock::now();

    time_end = chrono::high_resolution_clock::now();
    //cout<<"step 6 to 8 time: "<<(chrono::duration_cast<chrono::microseconds>(time_end - time_start)).count()<<endl;
    }



void my_fastbconv_m_tilde(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads)
{

    /*
    Require: Input in q
    Ensure: Output in Bsk U {m_tilde}
    */
 
    size_t base_Bsk_size = rns_tool->base_Bsk()->size(); // 14
    size_t base_q_size = rns_tool->base_q()->size();
    
    SEAL_ALLOCATE_GET_RNS_ITER(temp, (8192*4), base_q_size, pool); // Need to replace the hard coded value
    //multiply_poly_scalar_coeffmod(input, base_q_size, rns_tool->m_tilde().value(), rns_tool->base_q()->base(), temp);

    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < base_q_size; i++) {
        multiply_poly_scalar_coeffmod(input[i], temp.poly_modulus_degree(), rns_tool->m_tilde().value(), rns_tool->base_q()->base()[i], temp[i]);
    }
    // Now convert to Bsk
    //rns_tool->base_q_to_Bsk_conv()->fast_convert_array(temp, destination, pool);
    my_fast_convert_array(rns_tool->base_q_to_Bsk_conv(), temp, destination, pool, num_threads);

    // Finally convert to {m_tilde}
    //rns_tool->base_q_to_m_tilde_conv()->fast_convert_array(temp, destination + base_Bsk_size, pool);
    my_fast_convert_array(rns_tool->base_q_to_m_tilde_conv(), temp, destination + base_Bsk_size, pool, num_threads);
}


void my_fast_convert_array(const Pointer<BaseConverter> &conv, ConstRNSIter in, RNSIter out, MemoryPoolHandle pool, int num_threads)
{
    auto ibase_ = conv->ibase();
    auto obase_ = conv->obase();
    size_t ibase_size = ibase_.size();
    size_t obase_size = obase_.size();
    size_t count = in.poly_modulus_degree();

    // Note that the stride size is ibase_size
    SEAL_ALLOCATE_GET_STRIDE_ITER(temp, uint64_t, count, ibase_size, pool);
    
    omp_set_num_threads(num_threads);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < ibase_size;i++) {
        for(int j = 0; j < count;j++) {
            if (ibase_.inv_punctured_prod_mod_base_array()[i].operand == 1) {
                temp[j][i] = barrett_reduce_64(in[i][j], ibase_.base()[i]);
            } else {
                temp[j][i] = multiply_uint_mod(in[i][j], ibase_.inv_punctured_prod_mod_base_array()[i],ibase_.base()[i]);
            } 
           
        }
      
    }

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < obase_size; i++) {
        for(int j = 0; j < count; j++){
            out[i][j] = dot_product_mod(temp[j], conv->base_change_matrix()[i].get(), ibase_size, obase_.base()[i]);

        }       
    }

}

inline void my_multiply_poly_scalar_coeffmod(ConstRNSIter poly, std::size_t coeff_modulus_size, std::uint64_t scalar, ConstModulusIter modulus,
    RNSIter result, int num_threads)
{
    auto poly_modulus_degree = result.poly_modulus_degree();

    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < coeff_modulus_size; i++) {
        multiply_poly_scalar_coeffmod(poly[i], poly_modulus_degree, scalar, modulus[i], result[i]);
    }
}

void my_sm_mrq(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads)
{
    /*
    Require: Input in base Bsk U {m_tilde}
    Ensure: Output in base Bsk
    */

    size_t base_Bsk_size = rns_tool->base_Bsk()->size();;

    // The last component of the input is mod m_tilde
    ConstCoeffIter input_m_tilde = input[base_Bsk_size];
    const uint64_t m_tilde_div_2 = rns_tool->m_tilde().value() >> 1;

    auto coeff_count_ = rns_tool->coeff_count();
    // Compute r_m_tilde
    SEAL_ALLOCATE_GET_COEFF_ITER(r_m_tilde, coeff_count_, pool);
    multiply_poly_scalar_coeffmod( input_m_tilde, coeff_count_, rns_tool->neg_inv_prod_q_mod_m_tilde(), rns_tool->m_tilde(), r_m_tilde);

    omp_set_num_threads(num_threads);

    for(int i = 0; i < base_Bsk_size; i++) {
        MultiplyUIntModOperand prod_q_mod_Bsk_elt;
        prod_q_mod_Bsk_elt.set(rns_tool->prod_q_mod_Bsk()[i], rns_tool->base_Bsk()->base()[i]);
        SEAL_ITERATE(iter(input[i], r_m_tilde, destination[i]), coeff_count_, [&](auto J) {
            uint64_t temp = get<1>(J);
            if (temp >= m_tilde_div_2)
            {
                temp += rns_tool->base_Bsk()->base()[i].value() - rns_tool->m_tilde().value();
            }
            get<2>(J) = multiply_uint_mod(
                multiply_add_uint_mod(temp, prod_q_mod_Bsk_elt, get<0>(J),rns_tool->base_Bsk()->base()[i]), rns_tool->inv_m_tilde_mod_Bsk()[i],
                rns_tool->base_Bsk()->base()[i]);
        });
    }

 }

void my_fast_floor(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads)
{

    size_t base_q_size = rns_tool->base_q()->size();
    size_t base_Bsk_size = rns_tool->base_Bsk()->size();

    // Convert q -> Bsk
    my_fast_convert_array(rns_tool->base_q_to_Bsk_conv(),input, destination, pool, num_threads);

    // Move input pointer to past the base q components
    input += base_q_size;
    #pragma omp parallel for collapse(2) num_threads(num_threads)
    for(int i = 0; i < base_Bsk_size;i++) {
        for(int j = 0; j < rns_tool->coeff_count(); j++) {
            destination[i][j] = multiply_uint_mod(input[i][j] + (rns_tool->base_Bsk()->base()[i].value() - destination[i][j]), rns_tool->inv_prod_q_mod_Bsk()[i], rns_tool->base_Bsk()->base()[i]);
        }    
    }
}

void my_fastbconv_sk(const RNSTool *rns_tool, ConstRNSIter input, RNSIter destination, MemoryPoolHandle pool, int num_threads)
{
    /*
    Require: Input in base Bsk
    Ensure: Output in base q
    */

    size_t base_q_size = rns_tool->base_q()->size();
    size_t base_B_size = rns_tool->base_B()->size();

    auto coeff_count_ = rns_tool->coeff_count();
    // Fast convert B -> q; input is in Bsk but we only use B
    my_fast_convert_array(rns_tool->base_B_to_q_conv(), input, destination, pool, num_threads);


    // Compute alpha_sk
    // Fast convert B -> {m_sk}; input is in Bsk but we only use B
    SEAL_ALLOCATE_GET_COEFF_ITER(temp, coeff_count_, pool);
    my_fast_convert_array(rns_tool->base_B_to_m_sk_conv(), input, RNSIter(temp, coeff_count_), pool, num_threads);

    // Take the m_sk part of input, subtract from temp, and multiply by inv_prod_B_mod_m_sk_
    // Note: input_sk is allocated in input[base_B_size]
    SEAL_ALLOCATE_GET_COEFF_ITER(alpha_sk, coeff_count_, pool);

    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < coeff_count_;i++) {
         alpha_sk[i] = multiply_uint_mod(temp[i] + (rns_tool->m_sk().value() - input[base_B_size][i]), rns_tool->inv_prod_B_mod_m_sk(), rns_tool->m_sk());
    }

    // alpha_sk is now ready for the Shenoy-Kumaresan conversion; however, note that our
    // alpha_sk here is not a centered reduction, so we need to apply a correction below.
    const uint64_t m_sk_div_2 = rns_tool->m_sk().value() >> 1;


    #pragma omp parallel for num_threads(num_threads)
    for(int i = 0; i < base_q_size; i++) {

        MultiplyUIntModOperand prod_B_mod_q_elt;
        prod_B_mod_q_elt.set(rns_tool->prod_B_mod_q()[i], rns_tool->base_q()->base()[i]);

        MultiplyUIntModOperand neg_prod_B_mod_q_elt;
        neg_prod_B_mod_q_elt.set(rns_tool->base_q()->base()[i].value() - rns_tool->prod_B_mod_q()[i], rns_tool->base_q()->base()[i]);

        SEAL_ITERATE(iter(alpha_sk, destination[i]), coeff_count_, [&](auto J) {
            // Correcting alpha_sk since it represents a negative value
            if (get<0>(J) > m_sk_div_2)
            {
                get<1>(J) = multiply_add_uint_mod(
                    negate_uint_mod(get<0>(J), rns_tool->m_sk()), prod_B_mod_q_elt, get<1>(J), rns_tool->base_q()->base()[i]);
            }
            // No correction needed
            else
            {
                // It is not necessary for the negation to be reduced modulo the small prime
                get<1>(J) = multiply_add_uint_mod(get<0>(J), neg_prod_B_mod_q_elt, get<1>(J), rns_tool->base_q()->base()[i]);
            }
        });

    }

}


void my_relinearize_internal(SEALContext &context_, Ciphertext &encrypted, const RelinKeys &relin_keys, size_t destination_size, MemoryPoolHandle pool
, int num_threads) 
    {
        // Verify parameters.
        auto context_data_ptr = context_.get_context_data(encrypted.parms_id());
        if (!context_data_ptr)
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (relin_keys.parms_id() != context_.key_parms_id())
        {
            throw invalid_argument("relin_keys is not valid for encryption parameters");
        }

        size_t encrypted_size = encrypted.size();

        // Verify parameters.
        if (destination_size < 2 || destination_size > encrypted_size)
        {
            throw invalid_argument("destination_size must be at least 2 and less than or equal to current count");
        }
        if (relin_keys.size() < sub_safe(encrypted_size, size_t(2)))
        {
            throw invalid_argument("not enough relinearization keys");
        }

        // If encrypted is already at the desired level, return
        if (destination_size == encrypted_size)
        {
            return;
        }

        // Calculate number of relinearize_one_step calls needed
        size_t relins_needed = encrypted_size - destination_size;

        // Iterator pointing to the last component of encrypted
        auto encrypted_iter = iter(encrypted);
        encrypted_iter += encrypted_size - 1;

        SEAL_ITERATE(iter(size_t(0)), relins_needed, [&](auto I) {
            my_switch_key_inplace(context_,
                encrypted, *encrypted_iter, static_cast<const KSwitchKeys &>(relin_keys),
                RelinKeys::get_index(encrypted_size - 1 - I), pool, num_threads);
        });

        // Put the output of final relinearization into destination.
        // Prepare destination only at this point because we are resizing down
        encrypted.resize(context_, context_data_ptr->parms_id(), destination_size);
    }


    void my_switch_key_inplace( SEALContext &context_,
        Ciphertext &encrypted, ConstRNSIter target_iter, const KSwitchKeys &kswitch_keys, size_t kswitch_keys_index,
        MemoryPoolHandle pool, int num_threads)
    {
        auto parms_id = encrypted.parms_id();
        auto &context_data = *context_.get_context_data(parms_id);
        auto &parms = context_data.parms();
        auto &key_context_data = *context_.key_context_data();
        auto &key_parms = key_context_data.parms();
        auto scheme = parms.scheme();

        omp_set_num_threads(num_threads);
        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context_) || !is_buffer_valid(encrypted))
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (!target_iter)
        {
            throw invalid_argument("target_iter");
        }
        if (!context_.using_keyswitching())
        {
            throw logic_error("keyswitching is not supported by the context");
        }

        // Don't validate all of kswitch_keys but just check the parms_id.
        if (kswitch_keys.parms_id() != context_.key_parms_id())
        {
            throw invalid_argument("parameter mismatch");
        }

        if (kswitch_keys_index >= kswitch_keys.data().size())
        {
            throw out_of_range("kswitch_keys_index");
        }
        if (!pool)
        {
            throw invalid_argument("pool is uninitialized");
        }
        if (scheme == scheme_type::bfv && encrypted.is_ntt_form())
        {
            throw invalid_argument("BFV encrypted cannot be in NTT form");
        }
        if (scheme == scheme_type::ckks && !encrypted.is_ntt_form())
        {
            throw invalid_argument("CKKS encrypted must be in NTT form");
        }

        // Extract encryption parameters.
        size_t coeff_count = parms.poly_modulus_degree();
        size_t decomp_modulus_size = parms.coeff_modulus().size();
        auto &key_modulus = key_parms.coeff_modulus();
        size_t key_modulus_size = key_modulus.size();
        size_t rns_modulus_size = decomp_modulus_size + 1;
        auto key_ntt_tables = iter(key_context_data.small_ntt_tables());
        auto modswitch_factors = key_context_data.rns_tool()->inv_q_last_mod_q();

        // Size check
        if (!product_fits_in(coeff_count, rns_modulus_size, size_t(2)))
        {
            throw logic_error("invalid parameters");
        }

        // Prepare input
        auto &key_vector = kswitch_keys.data()[kswitch_keys_index];
        size_t key_component_count = key_vector[0].data().size();

        // Check only the used component in KSwitchKeys.
        for (auto &each_key : key_vector)
        {
            if (!is_metadata_valid_for(each_key, context_) || !is_buffer_valid(each_key))
            {
                throw invalid_argument("kswitch_keys is not valid for encryption parameters");
            }
        }

        // Create a copy of target_iter
        SEAL_ALLOCATE_GET_RNS_ITER(t_target, coeff_count, decomp_modulus_size, pool);
        set_uint(target_iter, decomp_modulus_size * coeff_count, t_target);

        // In CKKS t_target is in NTT form; switch back to normal form
        if (scheme == scheme_type::ckks)
        {
            inverse_ntt_negacyclic_harvey(t_target, decomp_modulus_size, key_ntt_tables);
        }

        // Temporary result
        auto t_poly_prod(allocate_zero_poly_array(key_component_count, coeff_count, rns_modulus_size, pool));

        #pragma omp parallel for
        for(int i = 0; i < rns_modulus_size; i++) {
            size_t key_index = (i == decomp_modulus_size ? key_modulus_size - 1 : i);
            // Product of two numbers is up to 60 + 60 = 120 bits, so we can sum up to 256 of them without reduction.
            size_t lazy_reduction_summand_bound = size_t(SEAL_MULTIPLY_ACCUMULATE_USER_MOD_MAX);
            size_t lazy_reduction_counter = lazy_reduction_summand_bound;

            // Allocate memory for a lazy accumulator (128-bit coefficients)
            auto t_poly_lazy(allocate_zero_poly_array(key_component_count, coeff_count, 2, pool));

            // Semantic misuse of PolyIter; this is really pointing to the data for a single RNS factor
            PolyIter accumulator_iter(t_poly_lazy.get(), 2, coeff_count);
            //#pragma omp parallel for num_threads(num_threads)
            for(int j = 0; j < decomp_modulus_size; j++) {
                SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);
                ConstCoeffIter t_operand;               
                
                // No need to perform RNS conversion (modular reduction)
                if (key_modulus[j] <= key_modulus[key_index])
                {
                    set_uint(t_target[j], coeff_count, t_ntt);
                }
                // Perform RNS conversion (modular reduction)
                else
                {
                    modulo_poly_coeffs(t_target[j], coeff_count, key_modulus[key_index], t_ntt);
                }
                // NTT conversion lazy outputs in [0, 4q)
                ntt_negacyclic_harvey_lazy(t_ntt, key_ntt_tables[key_index]);
                t_operand = t_ntt;
                
                // Multiply with keys and modular accumulate products in a lazy fashion
                SEAL_ITERATE(iter(key_vector[j].data(), accumulator_iter), key_component_count, [&](auto K) {
                    if (!lazy_reduction_counter)
                    {
                        SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)), coeff_count, [&](auto L) {
                            unsigned long long qword[2]{ 0, 0 };
                            multiply_uint64(get<0>(L), get<1>(L), qword);

                            // Accumulate product of t_operand and t_key_acc to t_poly_lazy and reduce
                            add_uint128(qword, get<2>(L).ptr(), qword);
                            get<2>(L)[0] = barrett_reduce_128(qword, key_modulus[key_index]);
                            get<2>(L)[1] = 0;
                        });
                    }
                    else
                    {
                        // Same as above but no reduction
                        SEAL_ITERATE(iter(t_operand, get<0>(K)[key_index], get<1>(K)), coeff_count, [&](auto L) {
                            unsigned long long qword[2]{ 0, 0 };
                            multiply_uint64(get<0>(L), get<1>(L), qword);
                            add_uint128(qword, get<2>(L).ptr(), qword);
                            get<2>(L)[0] = qword[0];
                            get<2>(L)[1] = qword[1];
                        });
                    }
                });

                if (!--lazy_reduction_counter)
                {
                    lazy_reduction_counter = lazy_reduction_summand_bound;
                }

            }

            // PolyIter pointing to the destination t_poly_prod, shifted to the appropriate modulus
            PolyIter t_poly_prod_iter(t_poly_prod.get() + (i * coeff_count), coeff_count, rns_modulus_size);

            // Final modular reduction
            //#pragma omp parallel for collapse(2)
            for(int k = 0; k < key_component_count; k++) {
                for(int l = 0; l < coeff_count; l++) {
                   if (lazy_reduction_counter == lazy_reduction_summand_bound)
                    {
                        (*t_poly_prod_iter[k])[l] = static_cast<uint64_t>(*accumulator_iter[k][l]);
                    }
                    else
                    {
                        (*t_poly_prod_iter[k])[l] = barrett_reduce_128(accumulator_iter[k][l].ptr(), key_modulus[key_index]);
                    }
                     
                }
            }

        }
        // Accumulated products are now stored in t_poly_prod

        // Perform modulus switching with scaling
        PolyIter t_poly_prod_iter(t_poly_prod.get(), coeff_count, rns_modulus_size);
        SEAL_ITERATE(iter(encrypted, t_poly_prod_iter), key_component_count, [&](auto I) {
            // Lazy reduction; this needs to be then reduced mod qi
            CoeffIter t_last(get<1>(I)[decomp_modulus_size]);
            inverse_ntt_negacyclic_harvey_lazy(t_last, key_ntt_tables[key_modulus_size - 1]);

            // Add (p-1)/2 to change from flooring to rounding.
            uint64_t qk = key_modulus[key_modulus_size - 1].value();
            uint64_t qk_half = qk >> 1;
            SEAL_ITERATE(t_last, coeff_count, [&](auto &J) {
                J = barrett_reduce_64(J + qk_half, key_modulus[key_modulus_size - 1]);
            });

            #pragma omp parallel for num_threads(num_threads)
            for(int j = 0; j < decomp_modulus_size; j++) {
                SEAL_ALLOCATE_GET_COEFF_ITER(t_ntt, coeff_count, pool);

                // (ct mod 4qk) mod qi
                uint64_t qi = key_modulus[j].value();
                if (qk > qi)
                {
                    // This cannot be spared. NTT only tolerates input that is less than 4*modulus (i.e. qk <=4*qi).
                    modulo_poly_coeffs(t_last, coeff_count, key_modulus[j], t_ntt);
                }
                else
                {
                    set_uint(t_last, coeff_count, t_ntt);
                }

                // Lazy substraction, results in [0, 2*qi), since fix is in [0, qi].
                uint64_t fix = qi - barrett_reduce_64(qk_half, key_modulus[j]);
                SEAL_ITERATE(t_ntt, coeff_count, [fix](auto &K) { K += fix; });

                uint64_t qi_lazy = qi << 1; // some multiples of qi
                inverse_ntt_negacyclic_harvey_lazy(get<1>(I)[j], key_ntt_tables[j]);


                // ((ct mod qi) - (ct mod qk)) mod qi
                SEAL_ITERATE(iter(get<1>(I)[j], t_ntt), coeff_count, [&](auto K) { get<0>(K) += qi_lazy - get<1>(K); });

                // qk^(-1) * ((ct mod qi) - (ct mod qk)) mod qi
                multiply_poly_scalar_coeffmod(get<1>(I)[j], coeff_count, modswitch_factors[j], key_modulus[j], get<1>(I)[j]);
                add_poly_coeffmod(get<1>(I)[j], get<0>(I)[j], coeff_count, key_modulus[j], get<0>(I)[j]);

            }

        });
    }

void my_bfv_multiply(SEALContext &context_, Ciphertext &encrypted1, Ciphertext &encrypted2, MemoryPoolHandle pool, int num_threads)
    {
        if (encrypted1.is_ntt_form() || encrypted2.is_ntt_form())
        {
            throw invalid_argument("encrypted1 or encrypted2 cannot be in NTT form");
        }

        // Extract encryption parameters.
        auto &context_data = *context_.get_context_data(encrypted1.parms_id());
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t base_q_size = parms.coeff_modulus().size();
        size_t encrypted1_size = encrypted1.size();
        size_t encrypted2_size = encrypted2.size();
        uint64_t plain_modulus = parms.plain_modulus().value();

        double new_scale = encrypted1.scale() * encrypted2.scale();

        // Check that scale is positive and not too large
        if (new_scale <= 0 || (static_cast<int>(log2(new_scale)) >= parms.plain_modulus().bit_count()))
        {
            throw invalid_argument("scale out of bounds");
        }
        omp_set_num_threads(num_threads);

        auto rns_tool = context_data.rns_tool();
        size_t base_Bsk_size = rns_tool->base_Bsk()->size();
        size_t base_Bsk_m_tilde_size = rns_tool->base_Bsk_m_tilde()->size();

        // Determine destination.size()
        size_t dest_size = sub_safe(add_safe(encrypted1_size, encrypted2_size), size_t(1));

        // Size check
        if (!product_fits_in(dest_size, coeff_count, base_Bsk_m_tilde_size))
        {
            throw logic_error("invalid parameters");
        }

        // Set up iterators for bases
        auto base_q = iter(parms.coeff_modulus());
        auto base_Bsk = iter(rns_tool->base_Bsk()->base());

        // Set up iterators for NTT tables
        auto base_q_ntt_tables = iter(context_data.small_ntt_tables());
        auto base_Bsk_ntt_tables = iter(rns_tool->base_Bsk_ntt_tables());


        // Microsoft SEAL uses BEHZ-style RNS multiplication. This process is somewhat complex and consists of the
        // following steps:
        //
        // (1) Lift encrypted1 and encrypted2 (initially in base q) to an extended base q U Bsk U {m_tilde}
        // (2) Remove extra multiples of q from the results with Montgomery reduction, switching base to q U Bsk
        // (3) Transform the data to NTT form
        // (4) Compute the ciphertext polynomial product using dyadic multiplication
        // (5) Transform the data back from NTT form
        // (6) Multiply the result by t (plain_modulus)
        // (7) Scale the result by q using a divide-and-floor algorithm, switching base to Bsk
        // (8) Use Shenoy-Kumaresan method to convert the result to base q

        // Resize encrypted1 to destination size
        encrypted1.resize(context_, context_data.parms_id(), dest_size);

        // This lambda function takes as input an IterTuple with three components:
        //
        // 1. (Const)RNSIter to read an input polynomial from
        // 2. RNSIter for the output in base q
        // 3. RNSIter for the output in base Bsk
        //
        // It performs steps (1)-(3) of the BEHZ multiplication (see above) on the given input polynomial (given as an
        // RNSIter or ConstRNSIter) and writes the results in base q and base Bsk to the given output
        // iterators.
        // auto behz_extend_base_convert_to_ntt = [&](auto I) {

            
        auto behz_extend_base_convert_to_ntt = [&](auto I) {

            SEAL_ALLOCATE_GET_RNS_ITER(temp, coeff_count, base_Bsk_m_tilde_size, pool);
            SEAL_ALLOCATE_GET_RNS_ITER(temp2, coeff_count, base_Bsk_m_tilde_size, pool);
            for(int j = 0; j < encrypted2_size; j++) {
                #pragma omp parallel for  num_threads(num_threads)  
                for(int i = 0; i < base_q_size ;i++) {
                    set_uint(get<0>(I)[i], coeff_count, get<1>(I)[i]);
                    ntt_negacyclic_harvey_lazy(get<1>(I)[i], base_q_ntt_tables[i]);
                    multiply_poly_scalar_coeffmod(get<0>(I)[i], temp2.poly_modulus_degree(), rns_tool->m_tilde().value(), rns_tool->base_q()->base()[i], temp2[i]);
                }
        
                my_fast_convert_array(rns_tool->base_q_to_Bsk_conv(), temp2, temp, pool, num_threads);
                my_fast_convert_array(rns_tool->base_q_to_m_tilde_conv(), temp2, temp + base_Bsk_size, pool, num_threads);

                ConstCoeffIter input_m_tilde = temp[base_Bsk_size];
                const uint64_t m_tilde_div_2 = rns_tool->m_tilde().value() >> 1;

                SEAL_ALLOCATE_GET_COEFF_ITER(r_m_tilde, rns_tool->coeff_count(), pool);
                multiply_poly_scalar_coeffmod( input_m_tilde, rns_tool->coeff_count(), rns_tool->neg_inv_prod_q_mod_m_tilde(), rns_tool->m_tilde(), r_m_tilde);

                #pragma omp parallel for num_threads(num_threads)   
                for(int i = 0; i < base_Bsk_size; i++) {
                    MultiplyUIntModOperand prod_q_mod_Bsk_elt;
                    prod_q_mod_Bsk_elt.set(rns_tool->prod_q_mod_Bsk()[i], rns_tool->base_Bsk()->base()[i]);
                    SEAL_ITERATE(iter(temp[i], r_m_tilde, get<2>(I)[i]), rns_tool->coeff_count(), [&](auto J) {
                        uint64_t temp = get<1>(J);
                        if (temp >= m_tilde_div_2)
                        {
                            temp += rns_tool->base_Bsk()->base()[i].value() - rns_tool->m_tilde().value();
                        }
                        get<2>(J) = multiply_uint_mod(
                            multiply_add_uint_mod(temp, prod_q_mod_Bsk_elt, get<0>(J),rns_tool->base_Bsk()->base()[i]), rns_tool->inv_m_tilde_mod_Bsk()[i],
                            rns_tool->base_Bsk()->base()[i]);
                    });
                    ntt_negacyclic_harvey_lazy(get<2>(I)[i], base_Bsk_ntt_tables[i]);

                }
            }
        };


        auto encrypted_iter2 = PolyIter(encrypted2);

        // Allocate space for a base q output of behz_extend_base_convert_to_ntt for encrypted1
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted1_q, encrypted1_size, coeff_count, base_q_size, pool);

        // Allocate space for a base Bsk output of behz_extend_base_convert_to_ntt for encrypted1
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted1_Bsk, encrypted1_size, coeff_count, base_Bsk_size, pool);

        // Perform BEHZ steps (1)-(3) for encrypted1
        SEAL_ITERATE(iter(encrypted1, encrypted1_q, encrypted1_Bsk), encrypted1_size, behz_extend_base_convert_to_ntt);



        // Repeat for encrypted2
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted2_q, encrypted2_size, coeff_count, base_q_size, pool);
        SEAL_ALLOCATE_GET_POLY_ITER(encrypted2_Bsk, encrypted2_size, coeff_count, base_Bsk_size, pool);
        SEAL_ITERATE(iter(encrypted2, encrypted2_q, encrypted2_Bsk), encrypted2_size, behz_extend_base_convert_to_ntt);


        //SEAL_ALLOCATE_GET_RNS_ITER(temp, coeff_count, base_Bsk_m_tilde_size, pool);
        //SEAL_ALLOCATE_GET_RNS_ITER(temp2, (8192*4), base_q_size, pool); // Need to replace the hard coded value

        // Allocate temporary space for the output of step (4)
        // We allocate space separately for the base q and the base Bsk components
        SEAL_ALLOCATE_ZERO_GET_POLY_ITER(temp_dest_q, dest_size, coeff_count, base_q_size, pool);
        SEAL_ALLOCATE_ZERO_GET_POLY_ITER(temp_dest_Bsk, dest_size, coeff_count, base_Bsk_size, pool);

        // Perform BEHZ step (4): dyadic multiplication on arbitrary size ciphertexts

        //#pragma omp parallel for
        for(int i = 0; i < dest_size;i++) {

            size_t curr_encrypted1_last = min<size_t>(i, encrypted1_size - 1);
            size_t curr_encrypted2_first = min<size_t>(i, encrypted2_size - 1);
            size_t curr_encrypted1_first = i - curr_encrypted2_first;
            // size_t curr_encrypted2_last = I - curr_encrypted1_last;

            // The total number of dyadic products is now easy to compute
            size_t steps = curr_encrypted1_last - curr_encrypted1_first + 1;

            // This lambda function computes the ciphertext product for BFV multiplication. Since we use the BEHZ
            // approach, the multiplication of individual polynomials is done using a dyadic product where the inputs
            // are already in NTT form. The arguments of the lambda function are expected to be as follows:
            //
            // 1. a ConstPolyIter pointing to the beginning of the first input ciphertext (in NTT form)
            // 2. a ConstPolyIter pointing to the beginning of the second input ciphertext (in NTT form)
            // 3. a ConstModulusIter pointing to an array of Modulus elements for the base
            // 4. the size of the base
            // 5. a PolyIter pointing to the beginning of the output ciphertext
            auto behz_ciphertext_product = [&](ConstPolyIter in1_iter, ConstPolyIter in2_iter,
                                                ConstModulusIter base_iter, size_t base_size, PolyIter out_iter) {
                // Create a shifted iterator for the first input
                auto shifted_in1_iter = in1_iter + curr_encrypted1_first;

                // Create a shifted reverse iterator for the second input
                auto shifted_reversed_in2_iter = reverse_iter(in2_iter + curr_encrypted2_first);

                // Create a shifted iterator for the output
                auto shifted_out_iter = out_iter[i];

                //#pragma omp parallel for collapse(2)
                for(int j = 0; j < steps; j++) {
                    #pragma omp parallel for num_threads(num_threads)
                    for(int k = 0; k < base_size; k++) {
                        SEAL_ALLOCATE_GET_COEFF_ITER(temp, coeff_count, pool);
                        dyadic_product_coeffmod(shifted_in1_iter[j][k], shifted_reversed_in2_iter[j][k], coeff_count, base_iter[k], temp);
                        add_poly_coeffmod(temp, shifted_out_iter[k], coeff_count, base_iter[k], shifted_out_iter[k]);

                    }
                }
            };

            // Perform the BEHZ ciphertext product both for base q and base Bsk
            behz_ciphertext_product(encrypted1_q, encrypted2_q, base_q, base_q_size, temp_dest_q);
            behz_ciphertext_product(encrypted1_Bsk, encrypted2_Bsk, base_Bsk, base_Bsk_size, temp_dest_Bsk);
        }


        // Perform BEHZ step (5): transform data from NTT form
        // Lazy reduction here. The following multiply_poly_scalar_coeffmod will correct the value back to [0, p)
        //inverse_ntt_negacyclic_harvey_lazy(temp_dest_q, dest_size, base_q_ntt_tables);
        //inverse_ntt_negacyclic_harvey_lazy(temp_dest_Bsk, dest_size, base_Bsk_ntt_tables);

        auto encrypted_iter1 = PolyIter(encrypted1);

        SEAL_ALLOCATE_GET_RNS_ITER(temp_q_Bsk, coeff_count, base_q_size + base_Bsk_size, pool);
        SEAL_ALLOCATE_GET_RNS_ITER(temp_Bsk, coeff_count, base_Bsk_size, pool);

        // Perform BEHZ steps (6)-(8)
        //#pragma omp parallel for 
        for(int i = 0; i < dest_size;i++) {
            #pragma omp parallel for num_threads(num_threads)
            for(int j = 0; j < temp_dest_Bsk.coeff_modulus_size();j++) { //14
                inverse_ntt_negacyclic_harvey_lazy(temp_dest_Bsk[i][j], base_Bsk_ntt_tables[j]);
                multiply_poly_scalar_coeffmod(temp_dest_Bsk[i][j], (temp_q_Bsk  + base_q_size).poly_modulus_degree(), plain_modulus, base_Bsk[j], (temp_q_Bsk + base_q_size)[j]);
                if(j < temp_dest_q.coeff_modulus_size()) {
                    inverse_ntt_negacyclic_harvey_lazy(temp_dest_q[i][j], base_q_ntt_tables[j]);
                    multiply_poly_scalar_coeffmod(temp_dest_q[i][j], temp_q_Bsk.poly_modulus_degree(), plain_modulus, base_q[j], temp_q_Bsk[j]);
                } 

            }
            my_fast_floor(rns_tool, temp_q_Bsk, temp_Bsk, pool, num_threads);
            my_fastbconv_sk(rns_tool, temp_Bsk, encrypted_iter1[i], pool, num_threads);

        }
        // Set the scale
    }


    void my_mod_switch_scale_to_next( SEALContext &context_,
        Ciphertext &encrypted, Ciphertext &destination, MemoryPoolHandle pool, int num_threads)
    {
        // Assuming at this point encrypted is already validated.
        auto context_data_ptr = context_.get_context_data(encrypted.parms_id());
        if (context_data_ptr->parms().scheme() == scheme_type::bfv && encrypted.is_ntt_form())
        {
            throw invalid_argument("BFV encrypted cannot be in NTT form");
        }
        if (context_data_ptr->parms().scheme() == scheme_type::ckks && !encrypted.is_ntt_form())
        {
            throw invalid_argument("CKKS encrypted must be in NTT form");
        }
        if (!pool)
        {
            throw invalid_argument("pool is uninitialized");
        }

        // Extract encryption parameters.
        auto &context_data = *context_data_ptr;
        auto &next_context_data = *context_data.next_context_data();
        auto &next_parms = next_context_data.parms();
        auto rns_tool = context_data.rns_tool();

        size_t encrypted_size = encrypted.size();
        size_t coeff_count = next_parms.poly_modulus_degree();
        size_t next_coeff_modulus_size = next_parms.coeff_modulus().size();

        Ciphertext encrypted_copy(pool);
        encrypted_copy = encrypted;
        auto encrypted_iter = iter(encrypted_copy);
        auto coeff_count_ = rns_tool->coeff_count();
        auto base_q_ = rns_tool->base_q();
        size_t base_q_size = base_q_->size();
        destination.resize(context_, next_context_data.parms_id(), encrypted_size);
        auto destination_iter = iter(destination);

        //#pragma omp parallel for   // Increases time, need to check again
        for(int i = 0; i < encrypted_size; i++) {
            //rns_tool->divide_and_round_q_last_inplace(encrypted_iter[i], pool);
            auto input = encrypted_iter[i];
            CoeffIter last_input = input[base_q_size - 1];

            // Add (qi-1)/2 to change from flooring to rounding
            Modulus last_modulus = (*base_q_)[base_q_size - 1];
            uint64_t half = last_modulus.value() >> 1;
            add_poly_scalar_coeffmod(last_input, coeff_count_, half, last_modulus, last_input);

            #pragma omp parallel for num_threads(num_threads)
            for(int j = 0; j < base_q_size - 1; j++) {
                SEAL_ALLOCATE_GET_COEFF_ITER(temp, coeff_count_, pool);
                modulo_poly_coeffs(last_input, coeff_count_, base_q_->base()[j], temp);
                uint64_t half_mod = barrett_reduce_64(half, base_q_->base()[j]);
                sub_poly_scalar_coeffmod(temp, coeff_count_, half_mod, base_q_->base()[j], temp);
                sub_poly_coeffmod(input[j], temp, coeff_count_, base_q_->base()[j], input[j]);
                multiply_poly_scalar_coeffmod(input[j], coeff_count_, rns_tool->inv_q_last_mod_q()[j], base_q_->base()[j], input[j]);                
            }

            set_poly(encrypted_iter[i], coeff_count, next_coeff_modulus_size, destination_iter[i]);

        }

        // Copy result to destination
        //destination.resize(context_, next_context_data.parms_id(), encrypted_size);
        // SEAL_ITERATE(iter(encrypted_copy, destination), encrypted_size, [&](auto I) {
        //     //set_poly(get<0>(I), coeff_count, next_coeff_modulus_size, get<1>(I));
        // });

        // Set other attributes
        destination.is_ntt_form() = encrypted.is_ntt_form();
    }


    void my_transform_to_ntt_inplace(SEALContext &context_, Ciphertext &encrypted, int num_threads) 
    {
        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context_) || !is_buffer_valid(encrypted))
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        auto context_data_ptr = context_.get_context_data(encrypted.parms_id());
        if (!context_data_ptr)
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (encrypted.is_ntt_form())
        {
            throw invalid_argument("encrypted is already in NTT form");
        }

        // Extract encryption parameters.
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t encrypted_size = encrypted.size();
        auto encrypted_iter = iter(encrypted);

        auto ntt_tables = iter(context_data.small_ntt_tables());

        // Size check
        if (!product_fits_in(coeff_count, coeff_modulus_size))
        {
            throw logic_error("invalid parameters");
        }

        #pragma omp parallel for collapse(2) num_threads(num_threads)
        for(int i = 0; i < encrypted_size; i++) {
            //#pragma omp parallel for
            for(int j = 0; j < encrypted.coeff_modulus_size(); j++) {
                ntt_negacyclic_harvey(encrypted_iter[i][j], ntt_tables[j]);
            }
        }

    

        // Finally change the is_ntt_transformed flag
        encrypted.is_ntt_form() = true;
    }


    void my_transform_from_ntt_inplace(SEALContext &context_, Ciphertext &encrypted_ntt, int num_threads) 
    {
        // Verify parameters.
        if (!is_metadata_valid_for(encrypted_ntt, context_) || !is_buffer_valid(encrypted_ntt))
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        auto context_data_ptr = context_.get_context_data(encrypted_ntt.parms_id());
        if (!context_data_ptr)
        {
            throw invalid_argument("encrypted_ntt is not valid for encryption parameters");
        }
        if (!encrypted_ntt.is_ntt_form())
        {
            throw invalid_argument("encrypted_ntt is not in NTT form");
        }

        // Extract encryption parameters.
        auto &context_data = *context_data_ptr;
        auto &parms = context_data.parms();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = parms.coeff_modulus().size();
        size_t encrypted_ntt_size = encrypted_ntt.size();
        auto encrypted_ntt_iter = iter(encrypted_ntt);

        auto ntt_tables = iter(context_data.small_ntt_tables());

        // Size check
        if (!product_fits_in(coeff_count, coeff_modulus_size))
        {
            throw logic_error("invalid parameters");
        }

        // Transform each polynomial from NTT domain
        //inverse_ntt_negacyclic_harvey(encrypted_ntt, encrypted_ntt_size, ntt_tables);

        #pragma omp parallel for collapse(2) num_threads(num_threads)
        for(int i = 0; i < encrypted_ntt_size; i++) {
            //#pragma omp parallel for
            for(int j = 0; j < encrypted_ntt.coeff_modulus_size(); j++) {
                inverse_ntt_negacyclic_harvey(encrypted_ntt_iter[i][j], ntt_tables[j]);
            }
        }

        // Finally change the is_ntt_transformed flag
        encrypted_ntt.is_ntt_form() = false;
    }

    void my_multiply_plain_ntt(SEALContext &context_, Ciphertext &encrypted_ntt, const Plaintext &plain_ntt, int num_threads)
    {
        // Verify parameters.
        if (!plain_ntt.is_ntt_form())
        {
            throw invalid_argument("plain_ntt is not in NTT form");
        }
        if (encrypted_ntt.parms_id() != plain_ntt.parms_id())
        {
            throw invalid_argument("encrypted_ntt and plain_ntt parameter mismatch");
        }

        // Extract encryption parameters.
        auto &context_data = *context_.get_context_data(encrypted_ntt.parms_id());
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t encrypted_ntt_size = encrypted_ntt.size();

        auto encrypted_ntt_iter = iter(encrypted_ntt);

        // Size check
        if (!product_fits_in(encrypted_ntt_size, coeff_count, coeff_modulus_size))
        {
            throw logic_error("invalid parameters");
        }

        double new_scale = encrypted_ntt.scale() * plain_ntt.scale();

        ConstRNSIter plain_ntt_iter(plain_ntt.data(), coeff_count);

        #pragma omp parallel for collapse(2) num_threads(num_threads)
        for( int i = 0; i < encrypted_ntt_size; i++) {
            for(int j = 0; j < coeff_modulus_size; j++) {
                dyadic_product_coeffmod(encrypted_ntt_iter[i][j], plain_ntt_iter[j], coeff_count, coeff_modulus[j], encrypted_ntt_iter[i][j]);
            }

        }

        // Set the scale
        encrypted_ntt.scale() = new_scale;
    }


    void my_rotate_internal(SEALContext context_, Ciphertext &encrypted, int steps, const GaloisKeys &galois_keys, MemoryPoolHandle pool, int num_threads)
    {
        auto context_data_ptr = context_.get_context_data(encrypted.parms_id());
        if (!context_data_ptr)
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (!context_data_ptr->qualifiers().using_batching)
        {
            throw logic_error("encryption parameters do not support batching");
        }
        if (galois_keys.parms_id() != context_.key_parms_id())
        {
            throw invalid_argument("galois_keys is not valid for encryption parameters");
        }

        // Is there anything to do?
        if (steps == 0)
        {
            return;
        }

        size_t coeff_count = context_data_ptr->parms().poly_modulus_degree();
        auto galois_tool = context_data_ptr->galois_tool();

        // Check if Galois key is generated or not.
        if (galois_keys.has_key(galois_tool->get_elt_from_step(steps)))
        {
            // Perform rotation and key switching
            my_apply_galois_inplace(context_, encrypted, galois_tool->get_elt_from_step(steps), galois_keys, move(pool), num_threads);
        }
        else
        {
            throw invalid_argument("Galois key not present");

        }
    }

    void my_conjugate_internal(SEALContext context_, Ciphertext &encrypted, const GaloisKeys &galois_keys, MemoryPoolHandle pool, int num_threads) 
    {
        // Verify parameters.
        auto context_data_ptr = context_.get_context_data(encrypted.parms_id());
        if (!context_data_ptr)
        {
            throw std::invalid_argument("encrypted is not valid for encryption parameters");
        }

        // Extract encryption parameters.
        auto &context_data = *context_data_ptr;
        if (!context_data.qualifiers().using_batching)
        {
            throw std::logic_error("encryption parameters do not support batching");
        }

        auto galois_tool = context_data.galois_tool();

        // Perform rotation and key switching
        my_apply_galois_inplace(context_, encrypted, galois_tool->get_elt_from_step(0), galois_keys, std::move(pool), num_threads);
    }


    void my_apply_galois_inplace(SEALContext context_, Ciphertext &encrypted, uint32_t galois_elt, const GaloisKeys &galois_keys, MemoryPoolHandle pool, int num_threads)
    {
        // Verify parameters.
        if (!is_metadata_valid_for(encrypted, context_) || !is_buffer_valid(encrypted))
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        // Don't validate all of galois_keys but just check the parms_id.
        if (galois_keys.parms_id() != context_.key_parms_id())
        {
            throw invalid_argument("galois_keys is not valid for encryption parameters");
        }

        auto &context_data = *context_.get_context_data(encrypted.parms_id());
        auto &parms = context_data.parms();
        auto &coeff_modulus = parms.coeff_modulus();
        size_t coeff_count = parms.poly_modulus_degree();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t encrypted_size = encrypted.size();
        // Use key_context_data where permutation tables exist since previous runs.
        auto galois_tool = context_.key_context_data()->galois_tool();

        // Size check
        if (!product_fits_in(coeff_count, coeff_modulus_size))
        {
            throw logic_error("invalid parameters");
        }

        // Check if Galois key is generated or not.
        if (!galois_keys.has_key(galois_elt))
        {
            throw invalid_argument("Galois key not present");
        }

        uint64_t m = mul_safe(static_cast<uint64_t>(coeff_count), uint64_t(2));

        // Verify parameters
        if (!(galois_elt & 1) || unsigned_geq(galois_elt, m))
        {
            throw invalid_argument("Galois element is not valid");
        }
        if (encrypted_size > 2)
        {
            throw invalid_argument("encrypted size must be 2");
        }

        if (parms.scheme() != scheme_type::bfv) {

            throw logic_error("scheme not implemented");
        }

        SEAL_ALLOCATE_GET_RNS_ITER(temp, coeff_count, coeff_modulus_size, pool);

        // DO NOT CHANGE EXECUTION ORDER OF FOLLOWING SECTION
        // BEGIN: Apply Galois for each ciphertext
        // Execution order is sensitive, since apply_galois is not inplace!

        auto encrypted_iter = iter(encrypted);
    
        #pragma omp parallel for num_threads(num_threads)
        for(int i = 0; i < coeff_modulus_size; i++) {
            galois_tool->apply_galois(encrypted_iter[0][i], galois_elt, coeff_modulus[i], temp[i]);
            set_poly(temp[i], coeff_count, 1, encrypted_iter[0][i]);
            galois_tool->apply_galois(encrypted_iter[1][i], galois_elt, coeff_modulus[i], temp[i]);
            set_zero_poly(coeff_count, 1, encrypted_iter[1][i]);

        }

        // Calculate (temp * galois_key[0], temp * galois_key[1]) + (ct[0], 0)
        my_switch_key_inplace(context_, encrypted, temp, static_cast<const KSwitchKeys &>(galois_keys), GaloisKeys::get_index(galois_elt), pool, num_threads);

    }
