
/* $Id$ */

#ifndef NCRESTARTDATA_H
#define NCRESTARTDATA_H

/* This header file will be used both in C++ and Fortran for preprocess */

/* dimensions */

#define _One_ "one"
#define _Six_ "six"
#define _Seven_ "seven"
#define _NAtoms_ "n_atoms"
#define _NAtoms3_ "n_atoms_3"
#define _NUse_ "n_use"
#define _MOLength_ "mo_length"
#define _ESPChargesLength_ "esp_charges_length"
#define _QMHessianLength_ "qm_hessian_length"

#define _SQPMNAtoms_ "sqpm_n_atoms"
#define _SQPMNImages_ "sqpm_n_images"
#define _SQPMNImagesPlus2_ "sqpm_n_images_2"
#define _SQPMNImages2NDim_ "sqpm_n_images_2_n_dim"
#define _SQPMNImagesNDim2_ "sqpm_n_images_n_dim_2"
#define _SQPM12NAtoms2_ "sqpm_12_natoms_2"
#define _SQPMNImages23NAtoms_ "sqpm_n_images_2_3_n_atoms"
#define _SQPMNImagesNAtoms_ "sqpm_n_images_n_atoms"

/* variables */

#define _NCRestartVersion_ "nc_restart_version"

#define _QMMMInformation_ "qmmm_information"
#define _QMMMGroup_ "qmmm_group"
#define _BoxDimensions_ "box_dimensions"
#define _AtomicNumbers_ "atomic_numbers"

#define _X_ "x"
#define _Y_ "y"
#define _Z_ "z"

#define _ActiveAtomList_ "active_atom_list"
#define _MOCoefficients_ "mo_coefficients"
#define _ESPCharges_ "esp_charges"
#define _QMHessian_ "qm_hessian"

#define _V_ "v"
#define _A_ "a"
#define _AOld_ "a_old"
#define _TimeStep_ "time_step"
#define _CurrentTimeStep_ "current_time_step"

#define _SQPMPath_ "sqpm_path"
#define _SQPMNewPath_ "sqpm_new_path"
#define _SQPMLastEnergies_ "sqpm_last_energies"
#define _SQPMLastGrads_ "sqpm_last_grads"
#define _SQPMTrustRadius_ "sqpm_trust_radius"
#define _SQPMLBPathH_ "sqpm_lb_path_H"
#define _SQPMPathFullXYZ_ "sqpm_path_full_xyz"
#define _SQPMRPXYZEnergyGradients "sqpm_rp_xyz_energy_gradients"
#define _SQPMPathXYZ_ "sqpm_path_xyz"
#define _SQPMPathESPCharges_ "sqpm_path_esp_charges"
#define _SQPMRPXYZEnergyGradients_ "sqpm_rp_xyz_energy_gradients"
#define _SQPMPathEnergyGradients_ "sqpm_path_energy_gradients"
#define _SQPMPathFullXYZ_ "sqpm_path_full_xyz"
#define _SQPMIterationStep_ "sqpm_iteration_step"

/* job type write */

#define _WriteQMMMInformation_                1
#define _WriteActiveAtomList_                 2
#define _WriteXYZ_                            3
#define _WriteMolecularOrbitalCoefficients_   4
#define _WriteSystemDimensions_               5
#define _WriteBoxDimensions_       _WriteSystemDimensions_
#define _WriteESPCharges_                     6
#define _WriteMolecularSystemInformation_     7
#define _WriteQMHessianFromOptimizeHess_      8
#define _WriteQMHessianFromFile_              9
#define _WriteMolecularDynamics_             10
#define _AddSQPMAtomsAndImagesNumber_        11
#define _WriteSQPM_                          12
#define _WriteSQPMLBPathH_                   13
#define _WriteSQPMRPXYZEnergyGradients_      14
#define _WriteSQPMPathFullXYZ_               15
#define _WriteSQPMPathESPCharges_            16
#define _WriteSQPMIterationStep_             17

/* job type read */

#define _CheckMolecularSystemInformation_   1001
#define _ReadActiveAtomList_                1002
#define _ReadXYZ_                           1003
#define _ReadMolecularDynamics_             1004
#define _ReadQMHessian_                     1005
#define _ReadMolecularOrbitalCoefficients_  1006
#define _ReadESPCharges_                    1007
#define _ReadSystemDimensions_              1008

#define _ReadSQPMAtomsAndImagesNumber_      1009
#define _ReadSQPM_                          1110
#define _ReadSQPMLBPathH_                   1111
#define _ReadSQPMRPXYZEnergyGradients_      1112
#define _ReadSQPMPathFullXYZ_               1113
#define _ReadSQPMPathESPCharges_            1114
#define _ReadSQPMIterationStep_             1115

/*  others */
#define _CreateNCRestart_                     0
#define _FlushNCRestart_                   9998
#define _CloseNCRestart_                   9999

#endif /* NCRESTART_H */
