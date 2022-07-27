#copyright (C) 2020  Ssohrab Borhanian
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import numpy as np
import gwbench.basic_relations as brs
from gwbench import network

############################################################################
### User Choices
############################################################################

# choose the desired detectors
network_spec = ['aLIGO_H', 'aLIGO_L', 'aLIGO_V']
#network_spec = ['ET_ET1', 'ET_ET2', 'ET_ET3']
#network_spec = ['CE2-40-CBO_C', 'CE2-40-CBO_N', 'CE2-40-CBO_S']

# initialize the network with the desired detectors
net = network.Network(network_spec)

# choose the desired waveform 
wf_model_name = 'heated_tf2_ecc'
# pass the chosen waveform to the network for initialization
net.set_wf_vars(wf_model_name=wf_model_name)

# pick the desired frequency range
#f = np.arange(5.,67.6,2**-4)

#define arrays for errors
err_Log_Mc = []
err_eta = []
err_H_eff5 = []
err_H_eff8 = []
err_ecc = []
snr_val = []
#chi_val = [] 
#ecc_val = []
corr_Log_Mc_ecc = []
corr_eta_ecc = []
corr_Heff5_ecc = []
corr_Heff8_ecc = []
corr_Log_Mc_eta = []
corr_Heff5_eta = []
corr_Heff8_eta = []
#input m1, m2 instead of Mc eta and then convert
#m1 = 15.
#m2 = 15.
chi = 0.5
ecc = 0.1
#Heff5_val = brs.Heff5(1, 1, m1, m2, chi, chi)
#Heff8_val = brs.Heff8(1, 1, m1, m2, chi, chi)
#print("Heff5 = {}, Heff8 = {}".format(Heff5_val, Heff8_val))

#loop for plotting
#j = 1

#while j>0:
#mass1 = [5.00025, 10.0005, 15.0008, 20.001, 25.0013, 30.0015, 35.0018, 40.002, 45.0023, 50.0025] #q=0.9999
#mass2 = [4.99975, 9.9995, 14.9992, 19.999, 24.9987, 29.9985, 34.9982, 39.998, 44.9977, 49.9975] #q=0.9999
#mass1 = [5., 10., 15., 20., 25., 30., 35., 40., 45., 50.] #q=1 
#mass2 = [5., 10., 15., 20., 25., 30., 35., 40., 45., 50.] #q=1
mass1 = [5.54772, 11.0954, 16.6432, 22.1909, 27.7386, 33.2863, 38.834, 44.3818, 49.9295, 55.4772] #q=1.24604
mass2 = [4.45228, 8.90456, 13.3568, 17.8091, 22.2614, 26.7137, 31.166, 35.6182, 40.0705, 44.5228] #q=1.24604
#mass1 = [6.66667, 13.3333, 20., 26.6667, 33.3333, 40., 46.6667, 53.3333, 60., 66.6667] #q=2
#mass2 = [3.33333, 6.66667, 10., 13.3333, 16.6667, 20., 23.3333, 26.6667, 30., 33.3333] #q=2
#mass1 = [7.5, 15., 22.5, 30., 37.5, 45., 52.5, 60., 67.5, 75.] #q=3
#mass2 = [2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25.]   #q=3
for m1, m2 in zip(mass1, mass2):
    Mc,eta = brs.Mc_eta_of_m1_m2(m1,m2)
    inj_params = {'Mc':    Mc,
                  'eta':   eta,
                  'chi1z': chi,
                  'chi2z': chi,
                  'DL':    200,
                  'tc':    0,
                  'phic':  0,
                  'iota':  np.pi/4,
                  'ra':    np.pi/4,
                  'dec':   np.pi/4,
                  'psi':   np.pi/4,
                  'gmst0': 0,
                  'Heff5': 1, #Heff5_val
                  'Heff8': 15, #Heff8_val
                  'e0': ecc
                 }
        #calculating isco frequency

        #Mc = inj_params["Mc"]
        #eta = inj_params["eta"]

    M = brs.M_of_Mc_eta(Mc,eta)
    f_isco = brs.f_isco_Msolar(M)


        #check the desired frequency range
    f = np.arange(10.,f_isco,2**-4)
#    f = np.arange(10.,100,2**-4)


        #wf_other_var_dic = {'Heff5': 0, 'Heff8': 0}
        # assign with respect to which parameters to take derivatives
    deriv_symbs_string = 'Mc eta Heff5 Heff8 e0'

        # assign which parameters to convert to cos or log versions
    conv_cos = ('iota','dec')
    conv_log = ('Mc','DL')

        # choose whether to take Earth's rotation into account
    use_rot = 0

        # pass all these variables to the network
    net.set_net_vars(
            f=f, inj_params=inj_params,
            deriv_symbs_string=deriv_symbs_string,
            conv_cos=conv_cos, conv_log=conv_log,
            use_rot=use_rot
            )

############################################################################
### GW benchmarking
############################################################################

        # compute the WF polarizations
    net.calc_wf_polarizations()
        # compute the WF polarizations and their derivatives
    net.calc_wf_polarizations_derivs_num()

        # setup antenna patterns, location phase factors, and PSDs
    net.setup_ant_pat_lpf_psds()

        # compute the detector responses
    net.calc_det_responses()
        # compute the detector responses and their derivatives
    net.calc_det_responses_derivs_num()

        # calculate the network and detector SNRs
    net.calc_snrs_det_responses()

        # calculate the network and detector Fisher matrices, condition numbers,
        # covariance matrices, error estimates, and inversion errors
    net.calc_errors()

        # calculate the 90%-credible sky area (in deg)
    net.calc_sky_area_90()

############################################################################
### Print results
############################################################################

        # print the contents of the detector objects (inside the network)
        #net.print_detectors()

        # print the contents of the network objects
    net.print_network()
        #print(net.errs)
    err_Log_Mc.append(net.errs["log_Mc"])
    err_eta.append(net.errs["eta"])
    err_H_eff5.append(net.errs["Heff5"])
    err_H_eff8.append(net.errs["Heff8"])
        #ecc_val.append(ecc)
        #chi_val.append(chi)
    err_ecc.append(net.errs["e0"])
    snr_val.append(net.snr)
    corr_Log_Mc_ecc.append(net.cov[0][4]/(net.errs["log_Mc"]*net.errs["e0"]))
    corr_eta_ecc.append(net.cov[1][4]/(net.errs["eta"]*net.errs["e0"]))
    corr_Heff5_ecc.append(net.cov[2][4]/(net.errs["Heff5"]*net.errs["e0"]))
    corr_Heff8_ecc.append(net.cov[3][4]/(net.errs["Heff8"]*net.errs["e0"]))
    corr_Log_Mc_eta.append(net.cov[0][1]/(net.errs["log_Mc"]*net.errs["eta"]))
    corr_Heff5_eta.append(net.cov[2][1]/(net.errs["Heff5"]*net.errs["eta"]))
    corr_Heff8_eta.append(net.cov[3][1]/(net.errs["Heff8"]*net.errs["eta"]))

        #print injection values of masses and spins
        #chi1,chi2 = inj_params["chi1z"],inj_params["chi2z"]
        #m1, m2 = brs.m1_m2_of_Mc_eta(Mc,eta)
        #print("m1 = {}, m2 = {}, chi1 = {}, chi2 = {} ".format(m1, m2, chi1, chi2))

        #ecc += 0.1
        #if ecc > 1: break
        #chi += 0.1
        #if chi>=0.95: break


#export outputs to a file
#----------------------------

with open ("/home/srishti.tiwari/gwbench/example_scripts/TidalHeatingEccResults/LIGO/error_q_1.24604_chi_0.5_e0_0.1.txt", "w") as f:
    f.write("#SNR\tm1\tm2\terr_Log_Mc\terr_eta\terr_H_eff5\terr_H_eff8\terr_ecc\tcorr_Mc_ecc\tcorr_eta_ecc\tcorr_Heff5_ecc\tcorr_Heff8_ecc\tcorr_Log_Mc_eta\tcorr_Heff5_eta\tcorr_Heff8_eta\n")
    for i in range(len(mass1)):
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n".format(snr_val[i], mass1[i], mass2[i], err_Log_Mc[i], err_eta[i], err_H_eff5[i],
                err_H_eff8[i], err_ecc[i], corr_Log_Mc_ecc[i], corr_eta_ecc[i], corr_Heff5_ecc[i], corr_Heff8_ecc[i], corr_Log_Mc_eta[i], corr_Heff5_eta[i], corr_Heff8_eta[i]))

#with open ("/home/srishti.tiwari/gwbench/example_scripts/TidalHeatingEccResults/error_m1_{}_m2_{}_chi_{}_e0_0.1to0.9.txt".format(m1, m2, chi), "w") as f:
        #f.write("ecc\terr_Log_Mc\terr_H_eff5\terr_H_eff8\terr_ecc\tcorr_Mc_ecc\tcorr_Heff5_ecc\tcorr_Heff8_ecc\n")
        #for i in range(0,len(ecc_val)):
                #f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(ecc_val[i], err_Log_Mc[i], err_H_eff5[i],
                         #err_H_eff8[i], err_ecc[i], corr_Log_Mc_ecc[i], corr_Heff5_ecc[i], corr_Heff8_ecc[i]))

