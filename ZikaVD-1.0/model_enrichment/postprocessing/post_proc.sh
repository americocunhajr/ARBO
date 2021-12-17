#!/bin/bash
rm *.dat *.m
cp ../outputData/sip_filtered_chain.m .
cp ../outputData/sip_raw_chain.m .
cp ../outputData/sip_filtered_chain_loglikelihood.m .
cp ../outputData/sip_raw_chain_logtarget.m .
#cp ../src/output50k/sip_filtered_chain.m .
#cp ../src/output200k/sip_filtered_chain.m .
#cp ../src/output500k/sip_filtered_chain.m .
#cp ../src/output1m/sip_filtered_chain.m .
#cp /workspace/repos/software/hydrogen/ip-catchall/src/outputData/sip_raw_chain.m .

cp ../outputData/sfp_qoi_seq.m .
#cp ../src/output50k/sfp_qoi_seq.m .
#cp ../src/output200k/sfp_qoi_seq.m .
#cp ../src/output500k/sfp_qoi_seq.m .
#cp ../src/output1m/sfp_qoi_seq.m .

./queso_m_to_dat

